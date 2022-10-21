"""
Cloudbutton Genomics Use case - Variant Caller Demo - processing fastq and
fasta inputs into chunks for map function to align all fastq x fasta combinations
using the gem3 mapper.
Reduce function merges mpileup and calls SNPs using SiNPle

Folder structure:
In bucket, create two folders:
1. fastqgz/  (for fastqgz file(s))
2. fasta/ (for fasta reference)
The program will then generate the following folders
3. fasta-chunks/ (for split fasta files)
4. fastq-indexes/ (for fastq index files)
5. outputs/ (pipeline output files - variant calling output from SiNPle)

USAGE:
command line examples

# 1. s3 commandline
# 1a. single-end test (to be updated)
python varcall_lithops_demo_v4.py -fq1 SRR6052133_1.fastq.gz -fa TriTrypDB-9.0_TbruceiTREU927_Genome_MbChr.fasta -b cloudbutton-variant-caller-input -nrfq 100000 -ncfa 270000 -itn 10 -rt lumimar/hutton-genomics-v04 -rtm 2048 > varcall_se.log
# 2a. paired end test (to be updated)
python varcall_lithops_demo_v4.py -fq1 SRR6052133_1.fastq.gz -fq2 SRR6052133_2.fastq.gz -fa TriTrypDB-9.0_TbruceiTREU927_Genome_MbChr.fasta -b cloudbutton-variant-caller-input -nrfq 400000 -ncfa 270000 -itn 10 -rt lumimar/hutton-genomics-v04 -rtm 2048 > varcall_pe.log
# 2. SRA commandline
python varcall_lithops_demo_v4.py -fq1 SRR6052133_1.fastq.gz -fq2 SRR6052133_2.fastq.gz -fa TriTrypDB-9.0_TbruceiTREU927_Genome_MbChr.fasta -b cloudbutton-variant-caller-input -nfq 100000 -nfa 270000 -itn 10 -ds SRA -fq SRR6052133 -rl 76 -rt lumimar/hutton-genomics-v03:12 -rtm 2048

TO DO:
Necessary
- correct coordinates for each chunk.
modify the iterdata to include the total number of chunks for each fastq file

Optional
- parallelise fastq.gz indexing (using a map-reduce iteration before the main map-reduce)


"""

###################################################################
###################################################################
# PACKAGES
# generic packages
import pathlib
import argparse
import math
import os.path
import time
import re
import sys
from random import randint

# specific packages (cloud computing, parallelisation etc)
import multiprocessing
import subprocess as sp
import boto3
import pandas as pd
from numpy import int64

# lithops functions and packages
from map_reduce_executor import MapReduce
from lithops.multiprocessing import util

# demo functions and packages
import lithopsgenetics
import lithopsgenetics.auxiliaryfunctions as af
import Metadata

###################################################################
###################################################################
##### PIPELINE SETTINGS
print("Variant Caller - Cloudbutton Genomics Use Case demo")
print("starting pipeline at "+ str(time.time()) + " - " + str(time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())))

# 0. DEBUGGING SETTINGS
gem_test=False
pre_processing_only=False
debug=True  #keep all af.printl outputs, if not False don't print them

print("DEBUG_TEST: running only gem indexer and mapper in map function: "+str(gem_test))
print("DEBUG_TEST: running only pre-processing: "+str(pre_processing_only))

# 1. PARSING COMMAND LINE
parser = argparse.ArgumentParser(description='Variant Caller - Cloudbutton Genomics Use Case demo')

# ADD COMMAND-LINE ARGUMENTS
# 1a. input file / names
# fastq in SRA database
parser.add_argument('-fq','--fq_seq_name', help='Fastq sequence name (for example SRR6052133) used for SRA database',required=False)

# fastq in s3 bucket
parser.add_argument('-fq1','--fastq1', help='Fastq file 1, stored in s3',required=False)
parser.add_argument('-fq2','--fastq2', help='Fastq file 2, stored in s3 (paired end sequencing) - optional',required=False)
parser.add_argument('-fa','--fasta',help='Fasta reference filename', required=True)

# 1b. input file location
parser.add_argument('-cl','--cloud_adr',help='cloud provider url prefix', required=False)
parser.add_argument('-b','--bucket',help='cloud provider bucket name', required=True)
parser.add_argument('-fb','--fbucket',help='cloud provider bucket name - for fasta file', required=False)

# 1c. fastq data source (SRA or s3)
parser.add_argument('-ds','--data_source',help='Data source', required=False)

# 1d. file splitting parameters
parser.add_argument('-nfq','--fastq_read_n', help='Number of reads per fastq chunk ',required=False)
# Modified, parser.add_argument('-nfa','--fasta_char_n',help='Number of characters per fasta chunk', required=True)
parser.add_argument('-nfa','--fasta_n_workers',help='Number of workers', required=True)
# Modified, parser.add_argument('-ofa','--fasta_char_overlap',help='bp overlap between fasta chunks', required=False)
parser.add_argument('-rl','--read_length',help='sequencing read length - used to calculate approx fastq chunk size', required=False)

# 1e. pipeline-specific parameters
parser.add_argument('-t','--tolerance',help='number of additional strata to include in filtration of map file', required=False)
parser.add_argument('-ff','--file_format',help='mpileup file format - csv or parquet', required=False)

# mpileup file format
# 1f. run settings
parser.add_argument('-itn','--iterdata_n',help='Number of iterdata elements to run', required=False)
parser.add_argument('-cf','--concur_fun',help='concurrent function quota limit', required=False)
parser.add_argument('-s3w','--temp_to_s3',help='Write intermediate temp files to s3 for debugging', required=False)
parser.add_argument('-rt','--runtime_id',help='runtime to use to execute the map-reduce', required=False)
parser.add_argument('-rtm','--runtime_mem',help='runtime memory to be assigned to each function - maximum 2048 MB', required=False)
parser.add_argument('-rtr','--runtime_memr',help='runtime memory to be assigned to reduce function - maximum 10240 MB', required=False)
parser.add_argument('-rts','--runtime_storage',help='runtime storage to be assigned to map function - maximum 10000 MB - currently set manually', required=False)
parser.add_argument('-bs','--buffer_size',help='memory in percentatge for buffer size - maximum 100%', required=False)
parser.add_argument('-ftm','--func_timeout_map',help='timeout for map function - maximum 900', required=False)
parser.add_argument('-ftr','--func_timeout_reduce',help='timeout for reduce function - maximum 900', required=False)
parser.add_argument('-sk','--skip_map',help='True/False; use mpileups generated by previous run, to run only reducer', required=False)
parser.add_argument('-lb','--loadbalancer',help='load balancer execution method: manual|select', required=False)

## PARSE COMMAND LINE ARGUMENTS
args = parser.parse_args()

# 2a. input file / names
fq_seqname = af.parse_optional_arg(args.fq_seq_name, None)
fastq_file = af.parse_optional_arg(args.fastq1, "")

# 2b. from input, determine whether it is paired- or single-end sequencing
fastq_file2, seq_type = af.parse_optional_arg(args.fastq2, "","paired-end", "single-end")
fasta_file = args.fasta

# 2c. input file location
cloud_adr = af.parse_optional_arg(args.cloud_adr, "aws")
BUCKET_NAME = args.bucket
FASTA_BUCKET = af.parse_optional_arg(args.fbucket, BUCKET_NAME)

# 2d. fastq data source (SRA)
datasource = af.parse_optional_arg(args.data_source, "s3")

# 3. file splitting parameters
# Fastq and fasta chunk sizes (fastq read no. multiplied by 4 to get number of lines)
fastq_read_n = int(af.parse_optional_arg(args.fastq_read_n, None))
fastq_chunk_size = 4*fastq_read_n  # used in the case of fastq stored in s3.
# Modified, fasta_chunk_size = int(args.fasta_char_n)
workers_fasta = int(args.fasta_n_workers)
# Modified, fasta_char_overlap = int(af.parse_optional_arg(args.fasta_char_overlap, 300))
read_length = int(af.parse_optional_arg(args.read_length, 150))

# 4. pipeline-specific parameters
tolerance = af.parse_optional_arg(args.tolerance, 0)
file_format = af.parse_optional_arg(args.file_format, "parquet")

# 5. run settings
iterdata_n, function_n = af.parse_optional_arg(args.iterdata_n, None,args.iterdata_n, "all")
concur_fun = int(af.parse_optional_arg(args.concur_fun, 10000))
temp_to_s3 = af.parse_optional_arg(args.temp_to_s3, False)
runtime_id = af.parse_optional_arg(args.runtime_id, 'lumimar/hutton-genomics-v03:18')
runtime_mem = af.parse_optional_arg(args.runtime_mem, 1024)
runtime_mem_r = af.parse_optional_arg(args.runtime_memr, 4096)
runtime_storage = af.parse_optional_arg(args.runtime_storage, 4000)
buffer_size = af.parse_optional_arg(args.buffer_size, "75%")
func_timeout_map = af.parse_optional_arg(args.func_timeout_map, 2400)
func_timeout_reduce = af.parse_optional_arg(args.func_timeout_reduce, 2400)
skip_map = af.parse_optional_arg(args.skip_map, False)
lb_method = af.parse_optional_arg(args.loadbalancer, "select")

# object prefixes (cloud folders)
fasta_folder = "fasta/"
fastq_folder = "fastqgz/"
fasta_folder_index = "fasta-indexs/"
idx_folder = "fastq-indexes/"
out_folder = "outputs/"
s3_temp_folder = "temp_outputs/"

# local log folders
local_log = "varcall_out/"
if not os.path.exists(local_log):
    os.makedirs(local_log)
lambda_log_stream = 'my_log_stream'


def mpileup_conversion(mpileup_file, fasta_chunk, fasta_key, fastq_chunk, file_format, BUCKET_NAME, storage):
        """
        Convert resulting data to csv/parquet and txt
        """
        
        # Filter mpileup file   
        with open(mpileup_file, 'r') as f:
            rows = f.read().splitlines()
            content = [row.split("\t") for row in rows]
            content.pop(-1) 
            del rows

        # Convert mpileup to Pandas dataframe
        df = pd.DataFrame(data=content)
        df.columns = df.columns.astype(str)
        df['1'] = df['1'].astype(int64)

        #create intermediate key
        fasta_key = fasta_chunks_prefix
        disallowed_characters = "._-!/·¡"
        for character in disallowed_characters:
            fasta_key = fasta_key.replace(character,"")

        # Create intermediate key
        fasta_chunk = fasta_chunk['id']+'_'+fasta_chunk['offset_base']
        max_index = df.iloc[-1]['1']
        intermediate_key = file_format + "/" + fasta_key + "_" + fasta_chunk + "-" + fastq_chunk[0] + "_chunk" + str(fastq_chunk[1]["number"]) + "_" + str(max_index)

        range_index = []
        x = 0
        while x < max_index:
            if(x < max_index):
                range_index.append(x)
            x = x + 100000

        x = x + (max_index - x)
        range_index.append(x)

        df3 = df.groupby(pd.cut(df['1'], range_index)).count()
        content = ""
        for i in range(len(df3)):
            content = content + str(range_index[i+1]) + ":" + str(df3.iloc[i,0]) + "\n"

        # Upload .txt file to storage
        storage.put_object(bucket=BUCKET_NAME, key=intermediate_key + ".txt", body=content)

        # Write the mpileup file to the tmp directory
        if file_format=="csv":
            df.to_csv(mpileup_file+".csv", index=False, header=False)
            with open(mpileup_file+".csv", 'rb') as f:
                storage.put_object(bucket=BUCKET_NAME, key=intermediate_key + ".csv", body=f)

        elif file_format=="parquet":
            df.to_parquet(mpileup_file+".parquet")
            with open(mpileup_file+".parquet", 'rb') as f:
                storage.put_object(bucket=BUCKET_NAME, key=intermediate_key + ".parquet", body=f)
        
        return [intermediate_key+"."+file_format, intermediate_key+".txt"]

def download_fastq(id, stage, fastq_chunk):
        """
        Download fastq chunks depending on source
        """   
        if seq_type == "paired-end":
            fastq1 = fastq_to_mapfun(fastq_chunk[0], fastq_chunk[1])
            base_name = os.path.splitext(fastq1)[0] #+'.pe'
            fastq2 = "yes"
        else:   # single-end sequencing
            fastq1 = fastq_to_mapfun(fastq_chunk[0], fastq_chunk[1])
            fastq2 = "no"
            base_name = os.path.splitext(fastq1)[0] #+'.se'
        return fastq1, fastq2, base_name


###################################################################
###################################################################
# MAP - REDUCE FUNCTIONS - MAP STAGE (STAGE A)
def map_alignment(id, fasta_chunk, fastq_chunk, storage):
    """
    map function to process fasta + fastq chunks to general sequence alignmments
    final format: mpileup
    """
    
    # CONTROL VARIABLES 
    stage="A"+str(randint(1000,9999))
    fasta_n=fasta_chunk['id']+'_'+str(fasta_chunk['offset_base'])
    fastq_n=re.sub('^[\s|\S]*number\':\s(\d*),[\s|\S]*$', r"\1", str(fastq_chunk))
    func_key=fq_seqname+"_fq"+fastq_n+"-"+fasta_n
    
    ###################################################################
    #### PROCESSING FASTQ CHUNKS
    ###################################################################
    # Download fastq chunk depending on source
    fastq1, fastq2, base_name = download_fastq(id, stage, fastq_chunk)

    ###################################################################
    #### PROCESSING FASTA CHUNKS
    ###################################################################
    fasta_folder_file = fasta_chunk['key_fasta'].split("/")
    fasta = af.copy_to_runtime(storage, FASTA_BUCKET, fasta_folder_file[0]+"/", fasta_folder_file[1], 
                                {'Range': f"bytes={fasta_chunk['offset_base']}-{fasta_chunk['last_byte+']}"}, fasta_chunk['offset_base'], fasta_chunk) # Download fasta chunk
    

    gem_ref_nosuffix = os.path.splitext(fasta)[0]
    gem_ref = gem_ref_nosuffix + '.gem'
    cpus=multiprocessing.cpu_count()
    sp.run(['gem-indexer', '--input', fasta, '--threads', str(cpus), '-o', gem_ref_nosuffix], capture_output=True)
    

    ###################################################################
    #### GENERATE ALIGNMENT AND ALIGNMENT INDEX (FASTQ TO MAP)
    ###################################################################
    sp.run(['/function/bin/map_index_and_filter_map_file_cmd_awsruntime.sh', gem_ref, fastq1, fastq2, base_name, datasource, seq_type], capture_output=True)

    # Reorganize file names
    map_index_file = base_name + "_map.index.txt"
    os.rename(map_index_file,"/tmp/" + func_key+ "_map.index.txt")
    map_index_file = "/tmp/" + func_key+ "_map.index.txt"
    old_filtered_map_file = base_name + "_filt_wline_no.map"
    filtered_map_file = base_name + "_" + str(id) + "_filt_wline_no.map"
    os.rename(old_filtered_map_file, filtered_map_file)
            
    # Copy intermediate files to storage for index correction
    map_index_file = af.copy_to_s3(storage, BUCKET_NAME, map_index_file, True, 'map_index_files/')
    filtered_map_file = af.copy_to_s3(storage, BUCKET_NAME, filtered_map_file, True, 'filtered_map_files/')
    map_index_file = map_index_file.replace("map_index_files/", "")
    filtered_map_file = filtered_map_file.replace("filtered_map_files/", "")

    return fasta_chunk, fastq_chunk, map_index_file, filtered_map_file, base_name, id


def map_alignment2(id, old_id, fasta_chunk, fastq_chunk, corrected_map_index_file, filtered_map_file, base_name, storage):
    """
    Second map  function, executed after the previous map function (map_alignment1) and the index correction.
    """
    
    ###################################################################
    #### RECOVER DATA FROM PREVIOUS MAP
    ###################################################################
    fasta_n=fasta_chunk['id']+'_'+str(fasta_chunk['offset_base'])
    fastq_n=re.sub('^[\s|\S]*number\':\s(\d*),[\s|\S]*$', r"\1", str(fastq_chunk))
    af.copy_to_runtime(storage, BUCKET_NAME, 'correctedIndex/', corrected_map_index_file)
    af.copy_to_runtime(storage, BUCKET_NAME, 'filtered_map_files/', filtered_map_file)
    corrected_map_index_file = "/tmp/" + corrected_map_index_file
    filtered_map_file = "/tmp/" + filtered_map_file
    # Modified, fasta_folder_file = fasta_chunk.split("/")
    fasta_folder_file = fasta_chunk['key_fasta'].split("/") 
    # Modified, fasta = af.copy_to_runtime(storage, FASTA_BUCKET, fasta_folder_file[0]+"/", fasta_folder_file[1], stage, id, debug) # Download fasta chunk
    fasta = af.copy_to_runtime(storage, FASTA_BUCKET, fasta_folder_file[0]+"/", fasta_folder_file[1], 
                                {'Range': f"bytes={fasta_chunk['offset_base']}-{fasta_chunk['last_byte+']}"}, fasta_chunk['offset_base'], fasta_chunk) # Download fasta chunk

    ###################################################################
    #### FILTER ALIGNMENTS (CORRECT .map FILE)
    ###################################################################
    map_filtering_output = sp.run(['/function/bin/map_file_index_correction_v3.sh', corrected_map_index_file, filtered_map_file, str(tolerance)], capture_output=True)  # change to _v3.sh and runtime 20

    corrected_map_file = base_name + "_" + str(old_id) + "_filt_wline_no_corrected.map"
        
    ###################################################################
    #### GENERATE MPILEUP FROM MAP FILE
    ###################################################################
    mpileup_out = sp.run(['/function/bin/gempileup_run.sh', corrected_map_file, fasta], capture_output=True)
    mpileup_file = corrected_map_file + ".mpileup"

        
    ###################################################################
    #### CONVERT MPILEUP TO PARQUET / CSV
    ###################################################################
    format_key, text_key = mpileup_conversion(mpileup_file, fasta_chunk, fasta_chunks_prefix, fastq_chunk, file_format, BUCKET_NAME, storage)

    return format_key, text_key 


def fastq_to_mapfun(fastq_file_key, fastq_chunk_data):
    '''
    Function executed within the map function to retrieve the relevant fastq chunk from object storage
    '''
    seq_name = fastq_file_key

    sp.call(['chmod','+x','fastq-dump'])
    sp.run(['vdb-config', '-i'])    # To supress a warning that appears the first time vdb-config is used

    # Report cloud identity so it can take data from s3 needed to be executed only once per vm
    output = str(sp.run(['vdb-config', '--report-cloud-identity', 'yes'], capture_output=True).stdout)
    
    os.chdir(f"/tmp")
    temp_fastq = f'/tmp/'+seq_name+f'_chunk{fastq_chunk_data["number"]}.fastq'
    data_output = sp.run(['fastq-dump', str(seq_name), '-X', str(int(fastq_chunk_data["start_line"])) , '-N', str(int(fastq_chunk_data["end_line"])),'-O',f'/tmp'],
                            capture_output=True)              
    os.rename(f'/tmp/'+seq_name+'.fastq', temp_fastq)
    
    return temp_fastq

###################################################################
###################################################################
# MAP - REDUCE FUNCTIONS - REDUCE STAGE 1 (STAGE B)

def create_sinple_key():

    fasta_chunk_n = math.ceil(int(iterdata_n) / len(fastq_list))
    keys = []
    fastq_file = fastq_list[0][0]

    i = 0
    for fasta_split in fasta_list:
        if i == fasta_chunk_n:
            break

        fasta_split = fasta_split.split("/")
        fasta_split = fasta_split[1]
        disallowed_characters = "._-!/·¡"
        for character in disallowed_characters:
            fasta_split = fasta_split.replace(character,"")

        fasta_split = fasta_split.split("split")
        fasta_file = fasta_split[0]
        fasta_split = fasta_split[1].split("fasta")

        fasta_file = fasta_file + "split_" + fasta_split[0] + ".fasta"
        # Modified, keys.append("multipart/" + fastq_file + "-" + fasta_file + "-" + seq_type + "_" +str(fasta_chunk_size) + "fa_" + str(fastq_chunk_size) + "fq_" + function_n + ".sinple")
        keys.append("multipart/" + fastq_file + "-" + fasta_file + "-" + seq_type + "_" +str(workers_fasta) + "fa_" + str(fastq_chunk_size) + "fq_" + function_n + ".sinple")
        i += 1

    fasta_file = fasta_file.split("split")

    # Modified, keys.append("multipart/" + fastq_file + "-" + fasta_file[0]+ ".fasta" + "-" + seq_type + "_" +str(fasta_chunk_size) + "fa_" + str(fastq_chunk_size) + "fq_" + function_n + ".sinple")
    keys.append("multipart/" + fastq_file + "-" + fasta_file[0]+ ".fasta" + "-" + seq_type + "_" +str(workers_fasta) + "fa_" + str(fastq_chunk_size) + "fq_" + function_n + ".sinple") 
    return keys


###################################################################
###################################################################
# MAP - REDUCE FUNCTIONS - REDUCE STAGE 2 (STAGE C)

def reduce_function(key, range, mpu_id, n_part, mpu_key, file_format, buffer_size, storage, id):
    stage="C"
    # check / prepare /tmp folder
    #af.file_and_folder_size("/tmp","at map function launch", stage, id, debug)
    # remove files from /tmp folder if they were left from previous runs of same script
    af.clear_tmp(stage, id, debug)
    if(os.path.exists('/tmp/reduce.mpileup')):
        os.remove('/tmp/reduce.mpileup')

    temp_mpileup = '/tmp/reduce.mpileup'
    fasta_n_first=re.sub('^\S*split_[0]*(\d*)\S*$', r"\1", key[0])
    af.printl("R\t"+fasta_n_first+"\t"+str(range['start']),stage,id,debug)

    tr0 = af.execution_time("retrieve mpileups","time_start",stage,id,debug)
    # Sequential
    # ---------------------------
    s3 = boto3.client('s3', endpoint_url="http://192.168.2.3:9000", aws_access_key_id="4OXIXRDDKRHB18U87QT2", aws_secret_access_key="gRphCMTppd9GS4CYIagGnwcttMZ0i3IOs46+O7XA")

    if file_format == "csv":
        expression = "SELECT * FROM s3object s WHERE cast(s._2 as int) BETWEEN %s AND %s" % (range['start'], range['end'])
        input_serialization = {'CSV': {}, 'CompressionType': 'NONE'}

    elif file_format == "parquet":
        expression = "SELECT * FROM s3object s WHERE s.\"1\" BETWEEN %s AND %s" % (range['start'], range['end'])
        input_serialization = {'Parquet': {}, 'CompressionType': 'NONE'}

    else:
        return "ERROR: Invalid format"
    
    for k in key:
        af.printl("s3 select key: "+str(k),stage,id,debug)
        fasta_n=re.sub('^\S*split_[0]*(\d*)\S*$', r"\1", k)
        fastq_n=re.sub('^\S*chunk(\d*)\S*$', r"\1", k)
        af.printl("R\t"+str(fasta_n)+"\t"+str(range['start'])+"\t"+fastq_n, stage,id,debug)
        #af.printl("fasta x fastq "+str(fasta_n)+"\t"+fastq_n, stage, id, debug)
        resp = s3.select_object_content(
            Bucket=BUCKET_NAME,
            Key=k,
            ExpressionType='SQL',
            Expression=expression,
            InputSerialization = input_serialization,
            OutputSerialization = {'CSV': {"FieldDelimiter" : "\t"}}
        )

        data = ""
        #record_count=0
        for event in resp['Payload']:
            if 'Records' in event:
                #record_count+=1
                records = event['Records']['Payload'].decode("UTF-8")
                #if record_count < 10:
                    #af.printl("record "+str(record_count) + " length: "+str(len(records)), stage, id, debug)
                data = data + records
        #af.printl("total number of mpileup records: "+str(record_count))

        #af.printl(data)

        wd = os.getcwd()
        os.chdir("/tmp")
        with open('/tmp/reduce.mpileup', 'a') as f:
            f.write(data)
        os.chdir(wd)
        del data
    # ---------------------------
    tr1 = af.execution_time("retrieve mpileups","time_end",stage,id,debug)
    af.printl(f' retrieve mpileups: execution_time:  {tr1 - tr0}: s',stage,id,debug)
    af.printl("r\t"+fasta_n_first+"\t"+str(range['start']),stage,id,debug)

    wd = os.getcwd()
    os.chdir("/tmp")
    
    tr2 = af.execution_time("merge mpileups","time_start",stage,id,debug)
    af.printl("S\t"+fasta_n_first+"\t"+str(range['start']),stage,id,debug)
    af.print_head_and_nlines(temp_mpileup, "mpileup file in reduce function", 10, stage, id, debug)
    # af.file_and_folder_size("/tmp","before mpileup_merge_reducev3.sh", stage, id, debug)
    af.copy_to_s3(stage, id, debug, storage, BUCKET_NAME, temp_mpileup, temp_to_s3, s3_temp_folder)
    af.printl("Starting Merge Script",stage,id,debug)
    sinple_out = sp.check_output(['bash', '/function/bin/mpileup_merge_reducev3_nosinple.sh', temp_mpileup, '/function/bin', buffer_size])
    sinple_out = sinple_out.decode('UTF-8')
    af.printl("Finished Merge Script",stage,id,debug)
    sinple_name=temp_mpileup+'_merged.mpileup'

    # write output to /tmp
    with open(sinple_name, 'w') as f:
        f.write(sinple_out)
    tr3 = af.execution_time("merge mpileups","time_end",stage,id,debug)
    af.printl("s\t"+fasta_n_first+"\t"+str(range['start']),stage,id,debug)
    af.printl(f' merge mpileups: execution_time:  {tr3 - tr2}: s',stage,id,debug)
    
    #af.file_and_folder_size("/tmp","after mpileup_merge_reducev3.sh", stage, id, debug)
    # copy output to s3
    af.copy_to_s3(stage, id, debug, storage, BUCKET_NAME, sinple_name, temp_to_s3, s3_temp_folder)

    os.chdir(wd)
    #Upload part
    s3 = boto3.client('s3', endpoint_url="http://192.168.2.3:9000", aws_access_key_id="4OXIXRDDKRHB18U87QT2", aws_secret_access_key="gRphCMTppd9GS4CYIagGnwcttMZ0i3IOs46+O7XA")
    part = s3.upload_part(
        Body = sinple_out,
        Bucket = BUCKET_NAME,
        Key = mpu_key,
        UploadId = mpu_id,
        PartNumber = n_part
    )
    tr4 = time.time()
    af.printl(f'reduce 2 - lambda: execution_time_total: {tr4 - tr0}: s',stage,id,debug)
    return {"PartNumber" : n_part, "ETag" : part["ETag"], "mpu_id": mpu_id}



if __name__ == "__main__":
    
    ### Preliminary steps:
    # 1. upload the fastq.gz file(s) into BUCKET_NAME/fastqgz
    # 2. Upload the fasta file into BUCKET_NAME/fasta
    
    run_id=str(randint(1000,9999))
    stage="PP"
    id="X"  # this is for function id, which is not present in this case
    PP_start = af.execution_time("preprocessing","time_start",stage,id,debug)

    ## RUN SETTINGS SUMMARY
    print("SERVERLESS VARIANT CALLER PIPELINE - USING LITHOPS + GEM-MAPPER + SiNPle")
    print("run id: "+run_id)
    
    print("RUN SETTINGS")
    print("command line: ")
    print(sys.argv[0]+str(sys.argv))
    print("\nBucket name: %s\n" % BUCKET_NAME )
    print("Sequencing type: " + seq_type + "\n")

    print("INPUT FILES")
    print("Fastq file 1: %s" % fastq_file )
    print("Fastq file 2: %s" % fastq_file2 )
    print("Fasta file: %s" % fasta_file )

    print("\nFILE SPLITTING SETTINGS")
    print("Fastq chunk size: %s lines" % str(fastq_chunk_size) )
    # Modified, print("Fasta chunk size: %s characters" % str(fasta_chunk_size) ) 
    print("Fasta number of workers: %s" % str(workers_fasta) ) 
    # Modified, print("Fasta overlap size: %s characters" % str(fasta_char_overlap) )

    print("\nOTHER RUN SETTINGS")
    if function_n == "all":
        print("Number of functions spawned: %s" % function_n )
    else:
        print("Number of functions spawned: %s" % iterdata_n )
    print("Runtime used: %s" % runtime_id )
    print("Runtime memory - map function: %s" % runtime_mem )
    print("Runtime memory - reduce function: %s" % runtime_mem_r )
    print("Runtime storage - map function: %s" % runtime_storage )
    print("Reduce load balancer method: %s" % lb_method )
    print("############\n")
    ###################################################################
    ###################################################################


    ###################################################################
    ###################################################################
    # 1. GENERATE LIST OF FASTQ CHUNKS (BYTE RANGES)
    t0 = af.execution_time("fastq list","time_start",stage,id,debug)

    num_spots = 0
    metadata = Metadata.SraMetadata()
    arr_seqs = [fq_seqname]
    if datasource == "SRA":
        accession = metadata.efetch_sra_from_accessions(arr_seqs)
        seq_type = accession['pairing'].to_string(index=False)
        num_spots = accession['spots'].to_string(index=False)
        fastq_size = accession['run_size'].to_string(index=False)
        print("Retrieving data from sra...")
        print("Sequence type: " + seq_type)
        print("Number of spots: " + num_spots)
        print("fastq size: " + fastq_size)

    # Generate fastq chunks
    fastq_list = lithopsgenetics.prepare_fastq(cloud_adr, BUCKET_NAME, idx_folder, fastq_folder, fastq_read_n, seq_type, fastq_file, fq_seqname, datasource, num_spots, fastq_file2)

    t1 = af.execution_time("fastq list","time_end",stage,id,debug)
    print(f'PP:0: fastq list: execution_time: {t1 - t0}: s')


    ###################################################################
    ###################################################################
    # 2. GENERATE LIST OF FASTA CHUNKS
    # It creates fasta chunks if they are not present
    t2 = af.execution_time("fasta list","time_start",stage,id,debug)

    # fasta file chunk prefix:
    # the prefix contains the name of the fasta reference minus the suffix,
    # followed by block length and then "_split_" i.e. for hg19.fa
    # with chunks with block length of 400000 it would be hg19_400000_split_
    # Modified,  fasta_chunks_prefix = re.sub("\.fasta|\.fa|\.fas", "_" + str(fasta_chunk_size) + "split_", fasta_file)
    fasta_chunks_prefix = pathlib.Path(fasta_file).stem
    # Generate fasta chunks
    # Modified,  fasta_list = lithopsgenetics.prepare_fasta(cloud_adr, runtime_id, FASTA_BUCKET, fasta_folder, fasta_file, fasta_folder_index, fasta_chunk_size, fasta_chunks_prefix, fasta_char_overlap)
    fasta_index = lithopsgenetics.prepare_fasta(FASTA_BUCKET, fasta_folder, fasta_file, fasta_folder_index, workers_fasta)
    
    t3 = af.execution_time("fasta list","time_end",stage,id,debug)
    print(f'PP:0: fasta list: execution_time: {t3 - t2}: s')


    ###################################################################
    ###################################################################
    # 3. GENERATE ITERDATA
    # iterdata will consist of a list with each fastq chunk paired with each fasta chunk
    print("\nGenerating iterdata")
    t4 = af.execution_time("iterdata","time_start",stage,id,debug)
    iterdata, fastq_set_n, num_chunks = lithopsgenetics.generate_alignment_iterdata(FASTA_BUCKET, fastq_list, fasta_index, fasta_folder+fasta_file, iterdata_n)

    # split iterdata if number higher than concurrent functions
    iterdata_n=len(iterdata)
            
    t5 = af.execution_time("iterdata","time_end",stage,id,debug)
    print(f'PP:0: iterdata: execution_time: {t5 - t4}: s')


    ###################################################################
    ###################################################################
    # 4. PREPROCESSING SUMMARY
    print("\nITERDATA LIST")
    print("number of fasta chunks: " + str(fastq_set_n))
    print("number of fastq chunks: " + str(len(fastq_list)))
    print("fasta x fastq chunks: "+ str(fastq_set_n*len(fastq_list)))
    print("number of iterdata elements: " + str(iterdata_n))
    PP_end = af.execution_time("preprocessing","time_start",stage,id,debug)
    print(f'PP:0: preprocessing: execution_time: {PP_end - PP_start}: s')

    if pre_processing_only==False:
        ###################################################################
        ###################################################################
        # 5. MAP-REDUCE
        stage="MR"
        id=0
        start = af.execution_time("map-reduce","time_start",stage,id,debug)
        mapreduce = MapReduce(map_alignment, map_alignment2, reduce_function, create_sinple_key, runtime_id, runtime_mem, runtime_mem_r, buffer_size, file_format, func_timeout_map, func_timeout_reduce, 'DEBUG', BUCKET_NAME, stage, id, debug, skip_map,10000,lb_method, fastq_set_n, run_id, iterdata_n, num_chunks, fq_seqname)
        map_time, creating_keys_time, reduce_time, multipart_upload_time = mapreduce(iterdata)
       
        print(str(map_time))
        exit()  # TODO

        
        end = af.execution_time("map-reduce","time_end",stage,id,debug)


        ###################################################################
        ###################################################################
        # 6. MAP-REDUCE SUMMARY
        
        print("MAP-REDUCE SUMMARY")
        af.printl("map phase: execution_time_total_varcall: "+str(map_time)+": s",stage, id, debug)
        af.printl("reduce 1 phase: execution_time_total_varcall: "+str(creating_keys_time)+": s",stage, id, debug)
        af.printl("reduce 2 phase: execution_time_total_varcall: "+str(reduce_time)+": s",stage, id, debug)
        af.printl("reduce 3 phase: execution_time_total_varcall: "+str(multipart_upload_time)+": s",stage, id, debug)
        af.printl("total time: execution_time_total_varcall: "+str(end - start)+": s",stage, id, debug)

        ###################################################################
        ###################################################################
