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
Preliminary steps:
    1. upload the fastq.gz file(s) into BUCKET_NAME/fastqgz
    2. Upload the fasta file into BUCKET_NAME/fasta

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
- modify the iterdata to include the total number of chunks for each fastq file

Optional
- parallelise fastq.gz indexing (using a map-reduce iteration before the main map-reduce)
"""

###################################################################
#### PACKAGES
###################################################################
# generic packages
import argparse
import os.path
import time
import re
import sys
from random import randint

# lithops functions and packages
from map_reduce_executor import MapReduce

# demo functions and packages
import lithopsgenetics
import lithopsgenetics.auxiliaryfunctions as af
import Metadata

# map/reduce functions
from map_reduce_functions import MapReduceFunctions

###################################################################
#### PIPELINE SETTINGS
###################################################################
parser = argparse.ArgumentParser(description='Variant Caller - Cloudbutton Genomics Use Case demo')

###################################################################
#### COMMAND-LINE ARGUMENTS
###################################################################
# File Names And Locations
# fastq in SRA database
parser.add_argument('-fq','--fq_seq_name', help='Fastq sequence name (for example SRR6052133) used for SRA database',required=False)
# fastq in s3 bucket
parser.add_argument('-fq1','--fastq1', help='Fastq file 1, stored in s3',required=False)
parser.add_argument('-fq2','--fastq2', help='Fastq file 2, stored in s3 (paired end sequencing) - optional',required=False)
parser.add_argument('-fa','--fasta',help='Fasta reference filename', required=True)
# input files locations
parser.add_argument('-cl','--cloud_adr',help='cloud provider url prefix', required=False)
parser.add_argument('-b','--bucket',help='cloud provider bucket name', required=True)
parser.add_argument('-fb','--fbucket',help='cloud provider bucket name - for fasta file', required=False)
# fastq data source (SRA or s3)
parser.add_argument('-ds','--data_source',help='Data source', required=False)

# File Splitting Parameters
parser.add_argument('-nfq','--fastq_read_n', help='Number of reads per fastq chunk ',required=False)
parser.add_argument('-nfa','--fasta_char_n',help='Number of characters per fasta chunk', required=True)
parser.add_argument('-ofa','--fasta_char_overlap',help='bp overlap between fasta chunks', required=False)

# Pipeline-Specific Parameters
parser.add_argument('-t','--tolerance',help='number of additional strata to include in filtration of map file', required=False)
parser.add_argument('-ff','--file_format',help='mpileup file format - csv or parquet', required=False)

# Run Settings
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

###################################################################
#### PARSE COMMAND LINE ARGUMENTS
###################################################################
args = parser.parse_args()

# FastQ names
fq_seqname = af.parse_optional_arg(args.fq_seq_name, None)
fastq_file = af.parse_optional_arg(args.fastq1, "")

# From input, determine whether it is paired- or single-end sequencing
fastq_file2, seq_type = af.parse_optional_arg(args.fastq2, "","paired-end", "single-end")
fasta_file = args.fasta

# Cloud Storage Settings
cloud_adr = af.parse_optional_arg(args.cloud_adr, "aws")
BUCKET_NAME = args.bucket
FASTA_BUCKET = af.parse_optional_arg(args.fbucket, BUCKET_NAME)

# Fastq data source (SRA)
datasource = af.parse_optional_arg(args.data_source, "s3")

# File Splitting Parameters
# Fastq and fasta chunk sizes (fastq read no. multiplied by 4 to get number of lines)
fastq_read_n = int(af.parse_optional_arg(args.fastq_read_n, None))
fastq_chunk_size = 4*fastq_read_n  # used in the case of fastq stored in s3.
fasta_chunk_size = int(args.fasta_char_n)
fasta_char_overlap = int(af.parse_optional_arg(args.fasta_char_overlap, 300))

# Pipeline-Specific Parameters
tolerance = af.parse_optional_arg(args.tolerance, 0)
file_format = af.parse_optional_arg(args.file_format, "parquet")

# Run Settings
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

# DEBUGGING SETTINGS
gem_test=False              # To execute only the gem indexer and mapper and skip the rest of map function
pre_processing_only=False   # To skip map/reduce
debug=True                  # Keep all af.printl outputs, if not False don't print them

print("DEBUG_TEST: running only gem indexer and mapper in map function: "+str(gem_test))
print("DEBUG_TEST: running only pre-processing: "+str(pre_processing_only))

# S3 prefixes (cloud folders)
fasta_folder = "fasta/"
fastq_folder = "fastqgz/"
split_fasta_folder = "fasta-chunks/"
idx_folder = "fastq-indexes/"
out_folder = "outputs/"
s3_temp_folder = "temp_outputs/"

# Local Log Folders
local_log = "varcall_out/"
if not os.path.exists(local_log):
    os.makedirs(local_log)
lambda_log_stream = 'my_log_stream'

 
if __name__ == "__main__":
    ###################################################################
    #### START THE PIPELINE
    ###################################################################
    run_id=str(randint(1000,9999))
    stage="PP"  
    id="X"      # this is for function id, which is not present in this case
    
    PP_start = af.execution_time("preprocessing","time_start",stage,id,debug)

    # Run settings summary
    print("Variant Caller - Cloudbutton Genomics Use Case demo")
    print("starting pipeline at "+ str(time.time()) + " - " + str(time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())))
    
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
    print("Fasta chunk size: %s characters" % str(fasta_chunk_size) )
    print("Fasta overlap size: %s characters" % str(fasta_char_overlap) )

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
    #### 1. GENERATE LIST OF FASTQ CHUNKS (BYTE RANGES)
    ###################################################################
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
    #### 2. GENERATE LIST OF FASTA CHUNKS (if not present)
    ###################################################################
    t2 = af.execution_time("fasta list","time_start",stage,id,debug)

    # Fasta File Chunk Prefix:
    # the prefix contains the name of the fasta reference minus the suffix,
    # followed by block length and then "_split_" i.e. for hg19.fa
    # with chunks with block length of 400000 it would be hg19_400000_split_
    fasta_chunks_prefix = re.sub("\.fasta|\.fa|\.fas", "_" + str(fasta_chunk_size) + "split_", fasta_file)
    
    # Generate fasta chunks
    fasta_list = lithopsgenetics.prepare_fasta(cloud_adr, runtime_id, FASTA_BUCKET, fasta_folder, fasta_file, split_fasta_folder, fasta_chunk_size, fasta_chunks_prefix, fasta_char_overlap)
    
    t3 = af.execution_time("fasta list","time_end",stage,id,debug)
    print(f'PP:0: fasta list: execution_time: {t3 - t2}: s')

    
    ###################################################################
    #### 3. GENERATE ITERDATA
    ###################################################################
    # The iterdata consists of an array where each element is a pair of a fastq chunk and a fasta chunk.
    # Since each fastq chunk needs to be paired with each fasta chunk, the total number of elements in
    # iterdata will be n_fastq_chunks * n_fasta_chunks.
    
    print("\nGenerating iterdata")
    t4 = af.execution_time("iterdata","time_start",stage,id,debug)
    
    # Generate iterdata
    iterdata, num_chunks = lithopsgenetics.generate_alignment_iterdata(fastq_list, fasta_list, iterdata_n)
    
    iterdata_n=len(iterdata)
    fastq_set_n=len(fasta_list) # number of fastq files aligned to each fasta chunk.
    
    # Split iterdata if number higher than concurrent functions
    iterdata_sets=[]
    iterdata_set=[]
    if iterdata_n > concur_fun:
        print("number of concurrent functions smaller than overall iterdata")
        count=0
        count_tot=0
        for el in iterdata:
            count+=1
            count_tot+=1
            fastq_count=el["fastq_chunk"][1]["number"]
            if ((count+fastq_set_n)>concur_fun and fastq_count==1):
                print("starting new iterdata set")
                print("count total: \t"+str(count_tot)+"\tcount: "+str(count)+"\tfastq_count "+str(fastq_count))
                print("length of finished set: "+ str(len(iterdata_set)))
                iterdata_sets.append(iterdata_set)
                iterdata_set=[]
                count=1
            iterdata_set.append(el)
            if count_tot==iterdata_n:
                iterdata_sets.append(iterdata_set)
        print("number of iterdata sets: "+str(len(iterdata_sets)))
        print("size of each iterdata set: ")
        for el in iterdata_sets:
            print(str(len(el)))
        
    t5 = af.execution_time("iterdata","time_end",stage,id,debug)
    print(f'PP:0: iterdata: execution_time: {t5 - t4}: s')

    
    ###################################################################
    #### 4. PREPROCESSING SUMMARY
    ###################################################################
    print("\nITERDATA LIST")
    print("number of fasta chunks: " + str(fastq_set_n))
    print("number of fastq chunks: " + str(len(fastq_list)))
    print("fasta x fastq chunks: "+ str(len(fasta_list)*len(fastq_list)))
    print("number of iterdata elements: " + str(iterdata_n))
    PP_end = af.execution_time("preprocessing","time_start",stage,id,debug)
    print(f'PP:0: preprocessing: execution_time: {PP_end - PP_start}: s')


    if pre_processing_only==False:
        ###################################################################
        #### 5. MAP-REDUCE
        ###################################################################
        stage="MR"
        id=0
        start = af.execution_time("map-reduce","time_start",stage,id,debug)
        functionobject = MapReduceFunctions(fq_seqname, datasource, seq_type, debug, BUCKET_NAME, fastq_folder, idx_folder, fasta_chunks_prefix, FASTA_BUCKET, gem_test, stage, file_format, tolerance)
        mapreduce = MapReduce(functionobject.map_alignment1, functionobject.map_alignment2, functionobject.reduce_function, functionobject.create_sinple_key, runtime_id, runtime_mem, runtime_mem_r, buffer_size, file_format, func_timeout_map, func_timeout_reduce, 'DEBUG', BUCKET_NAME, stage, id, debug, skip_map,10000,lb_method, fastq_set_n, run_id, iterdata_n, num_chunks, fq_seqname)
        map_time, creating_keys_time, reduce_time, multipart_upload_time = mapreduce(iterdata, iterdata_sets)
        end = af.execution_time("map-reduce","time_end",stage,id,debug)

        
        ###################################################################
        #### 6. MAP-REDUCE SUMMARY
        ###################################################################
        print("MAP-REDUCE SUMMARY")
        af.printl("map phase: execution_time_total_varcall: "+str(map_time)+": s",stage, id, debug)
        af.printl("reduce 1 phase: execution_time_total_varcall: "+str(creating_keys_time)+": s",stage, id, debug)
        af.printl("reduce 2 phase: execution_time_total_varcall: "+str(reduce_time)+": s",stage, id, debug)
        af.printl("reduce 3 phase: execution_time_total_varcall: "+str(multipart_upload_time)+": s",stage, id, debug)
        af.printl("total time: execution_time_total_varcall: "+str(end - start)+": s",stage, id, debug)
