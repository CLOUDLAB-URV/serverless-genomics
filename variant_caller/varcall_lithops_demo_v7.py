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
- redis function
- correct coordinates for each chunk.
modify the iterdata to include the total number of chunks for each fastq file

Optional
- parallelise fastq.gz indexing (using a map-reduce iteration before the main map-reduce)


"""
###################################################################
###################################################################
# PACKAGES
# generic packages
import argparse

import math
from math import floor
import os.path
import time
import re
import json
import sys
from pprint import pprint
from random import randint
# specific packages (cloud computing, parallelisation etc)
import redis
import multiprocessing
import subprocess as sp
import boto3
import shutil
import concurrent.futures
import pandas as pd
from numpy import int64
# lithops functions and packages
import lithops
from lithops import Storage
from map_reduce import MapReduce
import ec2_control.aws_ec2_control as ec2
from lithops.multiprocessing import util
# demo functions and packages
import lithopsgenetics
import lithopsgenetics.auxiliaryfunctions as af
import Metadata



###################################################################
###################################################################
# PIPELINE SETTINGS
print("Variant Caller - Cloudbutton Genomics Use Case demo")
print("starting pipeline at "+ str(time.time()) + " - " + str(time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())))

# 0. DEBUGGING SETTINGS
#threaded=True
gem_test=False
pre_processing_only=False
redis_index_correction=True
live_log=False
debug=True # keep all af.printl outputs, if not False don't print them
#print("DEBUG_TEST: running gem indexer with all available cpus: "+str(threaded))
print("DEBUG_TEST: running only gem indexer and mapper in map function: "+str(gem_test))
print("DEBUG_TEST: running only pre-processing: "+str(pre_processing_only))
print("DEBUG_TEST: redis index correction: "+str(redis_index_correction))
print("DEBUG_TEST: live log: "+str(redis_index_correction))
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
# 2. fastq data source (SRA or s3)
parser.add_argument('-ds','--data_source',help='Data source', required=False)
# 3. file splitting parameters
parser.add_argument('-nfq','--fastq_read_n', help='Number of reads per fastq chunk ',required=False)
parser.add_argument('-nfa','--fasta_char_n',help='Number of characters per fasta chunk', required=True)
parser.add_argument('-ofa','--fasta_char_overlap',help='bp overlap between fasta chunks', required=False)
parser.add_argument('-rl','--read_length',help='sequencing read length - used to calculate approx fastq chunk size', required=False)
# 4. pipeline-specific parameters
parser.add_argument('-t','--tolerance',help='number of additional strata to include in filtration of map file', required=False)
parser.add_argument('-ff','--file_format',help='mpileup file format - csv or parquet', required=False)
# mpileup file format
# 5. run settings
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
# timeout
# 6. Virtual Machine hosting Redis address
parser.add_argument('-ip','--ec2ip',help='IP address of the Virtual Machine hosting the Redis server', required=False)
parser.add_argument('-id','--ec2id',help='ID of the Virtual Machine hosting the Redis server', required=False)
parser.add_argument('-rg','--ec2region',help='Region the Virtual Machine hosting the Redis server is located', required=False)
parser.add_argument('-pem','--private_key',help='AWS private key', required=False)

# PARSE COMMAND LINE ARGUMENTS
args = parser.parse_args()
# 1a. input file / names
fq_seqname = af.parse_optional_arg(args.fq_seq_name, None)
fastq_file = af.parse_optional_arg(args.fastq1, "")
# From input, determine whether it is paired- or single-end sequencing
fastq_file2, seq_type = af.parse_optional_arg(args.fastq2, "","paired-end", "single-end")
fasta_file = args.fasta

# 1b. input file location
cloud_adr = af.parse_optional_arg(args.cloud_adr, "aws")
BUCKET_NAME = args.bucket
FASTA_BUCKET = af.parse_optional_arg(args.fbucket, BUCKET_NAME)

# 2. fastq data source (SRA)
datasource = af.parse_optional_arg(args.data_source, "s3")

# 3. file splitting parameters
# Fastq and fasta chunk sizes (fastq read no. multiplied by 4 to get number of lines)
fastq_read_n = int(af.parse_optional_arg(args.fastq_read_n, None))
fastq_chunk_size = 4*fastq_read_n  # used in the case of fastq stored in s3.
fasta_chunk_size = int(args.fasta_char_n)
fasta_char_overlap = int(af.parse_optional_arg(args.fasta_char_overlap, 300))
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
func_timeout_map = af.parse_optional_arg(args.func_timeout_map, 900)
func_timeout_reduce = af.parse_optional_arg(args.func_timeout_reduce, 900)
skip_map = af.parse_optional_arg(args.skip_map, False)
lb_method = af.parse_optional_arg(args.loadbalancer, "select")

# 6. Virtual Machine hosting Redis address
EC2_IP = af.parse_optional_arg(args.ec2ip, '54.146.89.181')
EC2_ID = af.parse_optional_arg(args.ec2id, 'i-0cee52f66655d990b')
REGION_NAME = af.parse_optional_arg(args.ec2region, 'us-east-1')
PEM = af.parse_optional_arg(args.private_key, 'aws-redis-damien.pem')
redis_password='lucio-redis-server'

# object prefixes (cloud folders)
fasta_folder = "fasta/"
fastq_folder = "fastqgz/"
split_fasta_folder = "fasta-chunks/"
idx_folder = "fastq-indexes/"
out_folder = "outputs/"
s3_temp_folder = "temp_outputs/"
# local folders
local_log = "varcall_out/"

# approximate fastq size
# calculation: size =(read_length*2.1) * number of reads
approx_fq_size_mb = floor((read_length * 2.1 * fastq_read_n)/100000)
approx_fa_size_mb = floor(int(fasta_chunk_size)/100000)
timestamp = time.strftime('%Y%m%d_%H%M', time.localtime())
if not os.path.exists(local_log):
    os.makedirs(local_log)
#sys.stdout = open(local_log+fq_seqname+"_"+str(approx_fq_size_mb)+"Mbfq_"+str(approx_fa_size_mb)+"Mbfa_"+str(iterdata_n)+"itn_"+timestamp+'.log','wt')
lambda_log_stream = 'my_log_stream' #
###################################################################
###################################################################

###################################################################
###################################################################
# MAP - REDUCE FUNCTIONS - MAP STAGE (STAGE A)
def map_alignment(id, fasta_chunk, fastq_chunk, storage):
    """
    map function to process fasta + fastq chunks to general sequence alignmments
    final format: mpileup
    """
    if (live_log==True):
        map_log = util.RemoteLogIOBuffer(lambda_log_stream)
        map_log.start()
    stage="A"+str(randint(1000,9999))
    af.printl("function number: "+ str(id),stage,id,debug)
    map_time_start = af.execution_time("map function","time_start",stage,id,debug)
    cpus=multiprocessing.cpu_count()
    af.printl("number of cpus: "+ str(cpus),stage,id,debug)
    # CONTROL VARIABLES 
    # fasta chunk ID

    fasta_n=re.sub(r'^\S*split_0*(\d*)\S*fasta', r'\1', fasta_chunk)
    fastq_n=re.sub('^[\s|\S]*number\':\s(\d*),[\s|\S]*$', r"\1", str(fastq_chunk))
    func_key=fq_seqname+"_fq"+fastq_n+"-fa"+fasta_n
    fq_set_key=fq_seqname+"_fq"+fastq_n
    af.printl("fasta chunk: "+str(fasta_chunk), stage, id, debug)
    af.printl("fasta number: "+str(fasta_n), stage, id, debug)
    af.printl("fastq number: "+str(fastq_n), stage, id, debug)
    af.printl("fastq split: "+str(func_key), stage, id, debug)
    
    # check / prepare /tmp folder
    #af.file_and_folder_size("/tmp","at map function launch", stage, id, debug)
    # remove files from /tmp folder if they were left from previous runs of same script
    af.clear_tmp(stage, id, debug)
    
    ###################################################################
    # PROCESSING FASTQ AND FASTA CHUNKS
    # 1. PROCESSING FASTQ CHUNK(S)
    ####start1
    print("E\t"+str(fasta_n)+"\t"+fastq_n)
    tm0 = af.execution_time("Fastq chunk","time_start",stage,id,debug)
    fastq1 =""
    fastq2 =""
    gem_ref =""
    af.printl("\n##########",stage,id,debug)
    if datasource == "s3":
        if seq_type == "paired-end":
            af.printl("processing paired-end fastq chunks",stage,id,debug)
            fastq1 = lithopsgenetics.fastq_to_mapfun("fastq1", fastq_chunk[0][0], fastq_chunk[0][1], BUCKET_NAME, fastq_folder, idx_folder, datasource, stage,id,debug)
            fastq2 = lithopsgenetics.fastq_to_mapfun("fastq2", fastq_chunk[1][0], fastq_chunk[1][1], BUCKET_NAME, fastq_folder, idx_folder, datasource, stage,id,debug)
            base_name = os.path.splitext(fastq1)[0]#+'.pe'
        else: # single-end sequencing
            af.printl("processing single-end fastq chunk",stage,id,debug)
            fastq1 = lithopsgenetics.fastq_to_mapfun("fastq", fastq_chunk[0], fastq_chunk[1], BUCKET_NAME, fastq_folder, idx_folder, datasource, stage,id,debug)
            base_name = os.path.splitext(fastq1)[0]+'.se'
            fastq2 = "no"
        # remove gz and gzi files (full-length fastq.gz)
        tmp_dir = os.listdir("/tmp")
        for item in tmp_dir:
            if item.endswith((".gz",".gzi")):
                os.remove(os.path.join(tmp_dir, item))
    elif datasource == "SRA":
        if seq_type == "paired-end":
            af.printl("processing paired-end fastq chunks",stage,id,debug)
            fastq1 = lithopsgenetics.fastq_to_mapfun("fastq", fastq_chunk[0], fastq_chunk[1], BUCKET_NAME, fastq_folder, idx_folder, datasource, stage,id,debug)
            base_name = os.path.splitext(fastq1)[0]#+'.pe'
            fastq2 = "yes"
        else: # single-end sequencing
            af.printl("processing single-end fastq chunk",stage,id,debug)
            fastq1 = lithopsgenetics.fastq_to_mapfun("fastq", fastq_chunk[0], fastq_chunk[1], BUCKET_NAME, fastq_folder, idx_folder, datasource, stage,id,debug)
            fastq2 = "no"
            base_name = os.path.splitext(fastq1)[0]#+'.se'
    
    af.print_head_and_nlines(fastq1, seq_type+" fastq file , "+ datasource, 8, stage, id, debug)
    if seq_type == "paired-end" and datasource == "s3":
        af.print_head_and_nlines(fastq2, seq_type+" fastq file 2 , "+ datasource, 8, stage, id, debug)
    tm1 = af.execution_time("Fastq chunk","time_end",stage,id,debug)
    af.printl(f' Fastq chunk: execution_time: {tm1 - tm0}: s',stage,id,debug)
    print("e\t"+str(fasta_n)+"\t"+fastq_n)
    ####end1

    # CHECKING IF GEM INDEX HAS ALREADY BEEN GENERATED
    gem_list = []
    gem_folder = "gem-chunks/"
    gem_ref=""
    gem_index_present=""
    #gem_chunks_prefix = re.sub("\.fasta|\.fa|\.fas", "_" + str(fasta_chunk_size) + "split_", fasta_file)
    try: 
        af.printl("Searching " + BUCKET_NAME + "/" + gem_folder + fasta_chunks_prefix, stage, id, debug)
        gem_list = storage.list_keys(BUCKET_NAME, prefix=gem_folder + fasta_chunks_prefix)
        af.printl("gem list: " + str(gem_list), stage, id, debug)
    except:
        af.printl("split gem folder empty / not found", stage, id, debug)

    # gem files has not been generated
    if not gem_list or gem_list == []:
        af.printl("processing fasta file and generating gem index in map function", stage, id, debug)

        # 2. PROCESSING FASTA CHUNK
        # 2a. copying fasta chunk to runtime
        ####start2a
        print("F\t"+str(fasta_n)+"\t"+fastq_n)
        tm2 = af.execution_time("Fasta chunk","time_start",stage,id,debug)
        fasta_chunk_folder_file = fasta_chunk.split("/")
        fasta = af.copy_to_runtime(storage, FASTA_BUCKET, fasta_chunk_folder_file[0]+"/", fasta_chunk_folder_file[1], stage, id, debug)
        af.print_head_and_nlines(fasta, "fasta file", 1, stage, id, debug)
        tm3 = af.execution_time("Fasta chunk","time_end",stage,id,debug)
        af.printl(f' Fasta chunk: execution_time: {tm3 - tm2}: s',stage,id,debug)
        ####end2a
        print("f\t"+str(fasta_n)+"\t"+fastq_n)
        
        af.file_and_folder_size("/tmp","fasta and fastq file sizes", stage, id, debug)

        # 2b. creating gem index file for fasta chunk [gem-indexer adds .gem to the output name]
        ####start2b
        print("G\t"+str(fasta_n)+"\t"+fastq_n)
        tm4 = af.execution_time("gem indexer","time_start",stage,id,debug)
        gem_ref_nosuffix = os.path.splitext(fasta)[0]
        af.printl("gem_ref_nosuffix " + gem_ref_nosuffix,stage,id,debug)
        gem_ref = gem_ref_nosuffix + '.gem'

        wd = os.getcwd()
        os.chdir("/tmp")
        indexer_output = sp.run(['gem-indexer', '--input', fasta, '--threads', str(cpus), '-o', gem_ref_nosuffix], capture_output=True)
        # if threaded ==False:
        #     cpus=1
        #     af.printl("gem indexing with 1 thread", stage, id, debug)
        #     indexer_output = sp.run(['gem-indexer', '--input', fasta, '--threads', str(cpus), '-o', gem_ref_nosuffix], capture_output=True)
        # else:
        #     af.printl("gem indexing with "+str(cpus)+" threads", stage, id, debug)
        #     indexer_output = sp.run(['gem-indexer', '--input', fasta, '--threads', str(cpus), '-o', gem_ref_nosuffix], capture_output=True)
        # Check size output of /tmp after gem
        files_gem = str(sp.check_output(["ls", '-la']))
        out = str(sp.check_output(["df", '-k']))
        result = re.search(r"(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\%\s+\/tmp", out)
        size, used, available, used_pct = (int(g) for g in result.groups())  
        af.printl("Base name: " + base_name,stage,id,debug)
        af.printl( '\n size: ' + str(size) + ' used: ' + str(used) + ' available: ' + str(available) + ' used_pct: ' + str(used_pct), stage, id, debug)
        af.print_indexer_output("INDEXER OUTPUT", indexer_output, "short", stage, id, debug)
        os.chdir(wd) 


        tm5 = af.execution_time("gem indexer","time_end",stage,id,debug)
        af.printl(f' gem indexer: execution_time: {tm5 - tm4}: s',stage,id,debug)
        print("g\t"+str(fasta_n)+"\t"+fastq_n)
        ####end2b


        af.copy_to_s3(stage, id, debug, storage, BUCKET_NAME, fastq1, temp_to_s3, s3_temp_folder)
        af.copy_to_s3(stage, id, debug, storage, BUCKET_NAME, fasta, temp_to_s3, s3_temp_folder)

        af.file_and_folder_size("/tmp","after gem indexing", stage, id, debug)

    else:
        gem_index_present=True
        for gem_file in gem_list:
            af.printl("fasta n before gem regex: "+str(fasta_n))
            query=r'' + str(fasta_n) + '.gem$'
            if re.search(query, gem_file):
                af.printl("found gem file "+gem_file+" associated with fasta chunk no. "+str(fasta_n))
                af.printl("retrieving gem index from s3", stage, id, debug)
                gem_file_name=fasta_chunks_prefix+str(fasta_n)+".gem"
                gem_ref = af.copy_to_runtime(storage, BUCKET_NAME, gem_folder, gem_file_name, stage, id, debug)
            else:
                af.printl("unmatched gem file "+gem_file)

        af.file_and_folder_size("/tmp","after gem retrieval", stage, id, debug)

    ###################################################################
    # 3. GENERATE ALIGNMENT AND ALIGNMENT INDEX (FASTQ TO MAP - gem3-mapper)
    ####start3
    tm6 = af.execution_time("gem mapper","time_start",stage,id,debug)
    print("H\t"+str(fasta_n)+"\t"+fastq_n)
    # the command below processes both single and paired end data
    af.printl("Base name: " + base_name,stage,id,debug)
    #af.print_head_and_nlines('/function/bin/map_index_and_filter_map_file_cmd_awsruntime.sh', "gem mapper bash script", 20, stage, id, debug)
    af.printl("datasource: "+datasource,stage,id,debug)
    mapper_output = sp.run(['/function/bin/map_index_and_filter_map_file_cmd_awsruntime.sh', gem_ref, fastq1, fastq2, base_name, datasource, seq_type], capture_output=True)
    af.print_formatted("MAPPER OUTPUT", mapper_output, stage, id, debug)

    tm7 = af.execution_time("gem mapper","time_end",stage,id,debug)
    af.printl(f' gem mapper: execution_time: {tm7 - tm6}: s',stage,id,debug)
    print("h\t"+str(fasta_n)+"\t"+fastq_n)
    ####end3

    # remove fastq and gem files from /tmp folder (no longer needed)
    af.printl("Deleting fastq and gem files from /tmp folder",stage,id,debug)
    os.remove(fastq1)
    if datasource == "s3" and seq_type == "paired-end":
        os.remove(fastq2)
    os.remove(gem_ref)

    af.file_and_folder_size("/tmp","after mapping", stage, id, debug)

    if gem_test == True:
        af.printl("running only gem indexer and mapper - skipping rest of map function", stage, id, debug)
    else: 
        
        map_index_file = base_name + "_map.index.txt"
        os.rename(map_index_file,"/tmp/" + func_key+ "_map.index.txt")
        map_index_file = "/tmp/" + func_key+ "_map.index.txt"
        filtered_map_file = base_name + "_filt_wline_no.map"
        af.print_head_and_nlines(map_index_file, "Map index file", 5, stage, id, debug)
        af.print_head_and_nlines(filtered_map_file, "Filtered map file", 5, stage, id, debug)
        #tmp_var=temp_to_s3
        #temp_to_s3=True
        af.copy_to_s3(stage, id, debug, storage, BUCKET_NAME, map_index_file, temp_to_s3, s3_temp_folder)
        af.copy_to_s3(stage, id, debug, storage, BUCKET_NAME, filtered_map_file, temp_to_s3, s3_temp_folder)
        #temp_to_s3=tmp_var

        corrected_map_index_file=""
        if redis_index_correction==True:
            ##################################################################
            # 2 REDIS
            tm8 = af.execution_time("redis","time_start",stage,id,debug)
            print("I\t"+str(fasta_n)+"\t"+fastq_n)

            r = redis.Redis(host='54.146.89.181', port=6379, password='lucio-redis-server',socket_connect_timeout=1)
            # SEND INDEX TO REDIS

            index_name= run_id+"__"+map_index_file
            af.printl("sending "+map_index_file+ " to redis",stage,id,debug)
            with open(map_index_file, 'r') as f:
                r.set(index_name, f.read())
                # af.printl("f.read output",stage,id,debug)
                # af.printl(f.read(), stage, id, debug)
                # af.printl("r.list_keys output", stage, id, debug)
                # for keys in r.list_keys(''):
                #     af.printl(str(keys),stage,id,debug)


            # SCRIPTS FOR MERGING INDICES IN REDIS
            # filter_merged_index.sh
            # merge_gem_alignment_metrics.sh

            # RECEIVE CORRECTED INDEX FROM REDIS
            item = fq_set_key+'.txt'
            af.printl("looking for corrected index "+item,stage,id,debug)
            split_list = r.keys(item)
            #af.printl("list of items in redis - before:",stage,id,debug)
            #for item in split_list:
            #    af.printl(item,stage,id,debug)
            found = False
            while not found :
                if len(split_list) > 0:
                    if split_list[0].decode('utf-8') == item:
                        found = True
                #af.printl("while loop - corrected index not found yet",stage,id,debug)
                time.sleep(0.2)
                split_list = r.keys(item)
                #af.printl("while loop - list of items in redis:\n" + str(split_list), stage, id, debug)
                #af.printl("item:\n" + str(item),stage,id, debug)
            #af.execution_time("redis","index_found",stage, id, debug)
            #split_list = r.list_keys('')
            #af.printl("list of items in redis - after:",stage,id,debug)
            #for item in split_list:
            #    af.printl(item,stage,id,debug)

            # process and check corrected map index
            corrected_map_index_obj = r.get(item)
            corrected_map_index_str = corrected_map_index_obj.decode('UTF-8')
            # if corrected index string is empty, then use the redis input index file.
            if corrected_map_index_str != "":
                corrected_map_index_file = base_name + "_map.corrected_index.txt"
                af.printl("corrected map index filename: "+corrected_map_index_file,stage,id,debug)
                with open(corrected_map_index_file, "w") as f:
                    f.write(corrected_map_index_str)

                af.print_head_and_nlines(corrected_map_index_file, "Corrected map index file", 10, stage, id, debug)
                tm9 = af.execution_time("redis","time_end", stage, id, debug)
                af.printl(f' redis: execution_time: {tm9 - tm8}: s',stage,id,debug)
                print("i\t"+str(fasta_n)+"\t"+fastq_n)

                af.copy_to_s3(stage, id, debug, storage, BUCKET_NAME, corrected_map_index_file, temp_to_s3, s3_temp_folder)
                

            
            else:

                corrected_map_index_file= map_index_file

            #exit(0)
            #print redis log in function.
            #redis_log = r.get_object('',"redis_log.txt")
            #af.print_formatted_redis("REDIS LOG", redis_log, stage, id, debug)


        else:
            corrected_map_index_file= map_index_file

        ##################################################################
        #3. FILTER ALIGNMENTS (CORRECT .map FILE)
        #map index file to be replaced with corrected map
        #index file

        tm10 = af.execution_time("correct_map","time_start", stage, id, debug)
        print("L\t"+str(fasta_n)+"\t"+fastq_n)
        map_filtering_output = sp.run(['/function/bin/map_file_index_correction.sh', corrected_map_index_file, filtered_map_file, str(tolerance)], capture_output=True)  # change to _v3.sh and runtime 20
        af.print_formatted("MAP FILTERING OUTPUT", map_filtering_output, stage, id, debug)
        tm11 = af.execution_time("correct_map","time_end", stage, id, debug)
        af.printl(f' correct_map: execution_time: {tm11 - tm10}: s',stage,id,debug)
        print("l\t"+str(fasta_n)+"\t"+fastq_n)

        corrected_map_file = base_name + "_filt_wline_no_corrected.map"
        # remove original non-corrected map file
        os.remove(filtered_map_file)

        af.file_and_folder_size("/tmp","after map correction", stage, id, debug)
        af.print_head_and_nlines(corrected_map_file, "Corrected map file", 10, stage, id, debug)
        #af.copy_to_s3(stage, id, debug, storage, BUCKET_NAME, corrected_map_file, temp_to_s3, s3_temp_folder)
        af.copy_to_s3(stage, id, debug, storage, BUCKET_NAME, corrected_map_file, "True", s3_temp_folder)
        
        ###################################################################
        # in the case of pre-loaded gem index, load fasta file for mpileup processing
        if gem_index_present:
            af.printl("processing fasta file for mpileup generation", stage, id, debug)

            # 2. PROCESSING FASTA CHUNK
            # 2a. copying fasta chunk to runtime
            ####start2a
            print("F\t"+str(fasta_n)+"\t"+fastq_n)
            tm2 = af.execution_time("Fasta chunk","time_start",stage,id,debug)
            fasta_chunk_folder_file = fasta_chunk.split("/")
            fasta = af.copy_to_runtime(storage, FASTA_BUCKET, fasta_chunk_folder_file[0]+"/", fasta_chunk_folder_file[1], stage, id, debug)
            af.print_head_and_nlines(fasta, "fasta file", 1, stage, id, debug)
            tm3 = af.execution_time("Fasta chunk","time_end",stage,id,debug)
            af.printl(f' Fasta chunk: execution_time: {tm3 - tm2}: s',stage,id,debug)
            ####end2a
            print("f\t"+str(fasta_n)+"\t"+fastq_n)
            
            af.file_and_folder_size("/tmp","fasta file size", stage, id, debug)



        ###################################################################
        # 4. GENERATE MPILEUP FROM MAP FILE
        tm12 = af.execution_time("gempileup","time_start", stage, id, debug)
        print("M\t"+str(fasta_n)+"\t"+fastq_n)
        mpileup_out = sp.run(['/function/bin/gempileup_run.sh', corrected_map_file, fasta], capture_output=True)
        af.print_formatted("MPILEUP OUTPUT", mpileup_out, stage, id, debug)
        tm13 = af.execution_time("gempileup","time_end", stage, id, debug)
        print("m\t"+str(fasta_n)+"\t"+fastq_n)
        af.printl(f' gempileup: execution_time: {tm13 - tm12}: s',stage,id,debug)
        
        mpileup_file = corrected_map_file + ".mpileup"
        os.remove(corrected_map_file)
        os.remove(fasta)
        af.print_head_and_nlines(mpileup_file, "mpileup file", 10, stage, id, debug)
        af.file_and_folder_size("/tmp","after gempileup", stage, id, debug)

        ###################################################################
        # FIX MPILEUP COORDINATES (pending on new fasta file splitting implementation)

        ###################################################################
        # 6. CONVERT MPILEUP TO PARQUET / CSV
        tm14 = af.execution_time("csv/parquet","time_start", stage, id, debug)
        print("N\t"+str(fasta_n)+"\t"+fastq_n)
        with open(mpileup_file, 'r') as f:
            #rows = f.readlines()
            rows = f.read().splitlines()
            content = [row.split("\t") for row in rows]
            content.pop(-1)  # modify to check if last line is empty
            del rows

        #Convert mpileup to Pandas dataframe -> '1' index position -> int
        df = pd.DataFrame(data=content)
        df.columns = df.columns.astype(str)
        df['1'] = df['1'].astype(int64)

        #create intermediate key
        fasta_chunk = fasta_chunk.split("_")
        fasta_key = fasta_chunks_prefix
        disallowed_characters = "._-!/·¡"
        for character in disallowed_characters:
            fasta_key = fasta_key.replace(character,"")

        df2 = df.iloc[-1]
        max_index = df2['1']



        # intermediate_key = file_format + "/" + fasta_key + "_" + fasta_chunk[-1] + "-" + fastq_chunk[0] + "_chunk" + str(fastq_chunk[1]["number"]) + "_" + str(max_index) + "." + file_format
        # af.printl("intermediate key: "+intermediate_key, stage, id )
        intermediate_key = file_format + "/" + fasta_key + "_" + fasta_chunk[-1] + "-" + fastq_chunk[0] + "_chunk" + str(fastq_chunk[1]["number"]) + "_" + str(max_index)
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
        #-----------------

        af.printl("intermediate key: \n"+ file_format + "/" + fasta_key + "_" + fasta_chunk[-1] + "-" + fastq_chunk[0] + "_chunk" + str(fastq_chunk[1]["number"]), stage, id )

        #write the mpileup file to the tmp directory
        map_output_file =""
        if file_format=="csv":
            map_output_file="csv/"+mpileup_file+".csv"
            df.to_csv(mpileup_file+".csv", index=False, header=False)
            #Upload the file to the s3
            with open(mpileup_file+".csv", 'rb') as f:
                storage.put_object(bucket=BUCKET_NAME, key=intermediate_key + ".csv", body=f)

            storage.put_object(bucket=BUCKET_NAME, key=intermediate_key + ".txt", body=content)
        elif file_format=="parquet":
            map_output_file="parquet/"+mpileup_file+".parquet"
            df.to_parquet(mpileup_file+".parquet")
            #Upload the file to the s3
            with open(mpileup_file+".parquet", 'rb') as f:
                storage.put_object(bucket=BUCKET_NAME, key=intermediate_key + ".parquet", body=f)
            storage.put_object(bucket=BUCKET_NAME, key=intermediate_key + ".txt", body=content)
        else:
            af.printl("file format not supported: "+file_format,stage,id,debug)

        tm15 = af.execution_time("csv/parquet","time_end", stage, id, debug)
        print("n\t"+str(fasta_n)+"\t"+fastq_n)
        af.printl(f' csv/parquet: execution_time: {tm15 - tm14}: s',stage,id,debug)
        af.printl(f' map - lambda: execution_time_total: {tm15 - map_time_start}: s',stage,id,debug)
        af.file_and_folder_size("/tmp","after conversion to csv/parquet", stage, id, debug)
        if (live_log==True):
            map_log.stop()
        #return intermediate_key
        return [intermediate_key+"."+file_format, intermediate_key+".txt"]


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
        keys.append("multipart/" + fastq_file + "-" + fasta_file + "-" + seq_type + "_" +str(fasta_chunk_size) + "fa_" + str(fastq_chunk_size) + "fq_" + function_n + ".sinple")
        i += 1

    fasta_file = fasta_file.split("split")

    keys.append("multipart/" + fastq_file + "-" + fasta_file[0]+ ".fasta" + "-" + seq_type + "_" +str(fasta_chunk_size) + "fa_" + str(fastq_chunk_size) + "fq_" + function_n + ".sinple")
    return keys


###################################################################
###################################################################
# MAP - REDUCE FUNCTIONS - REDUCE STAGE 2 (STAGE C)

def reduce_function(key, range, mpu_id, n_part, mpu_key, file_format, buffer_size, storage, id):
    if (live_log==True):
        red1_log = util.RemoteLogIOBuffer(lambda_log_stream)
        red1_log.start()
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
    s3 = boto3.client('s3')

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

    chr_table = af.copy_to_runtime(storage, FASTA_BUCKET, "cache/", "hashtable.data", stage, id, debug)
    af.printl("chromosome hash table: "+str(chr_table), stage, id, debug)
    with open(chr_table, 'r') as f:
        for line in f:
            af.printl(line, stage, id, debug)

    wd = os.getcwd()
    os.chdir("/tmp")
    
    tr2 = af.execution_time("merge mpileups","time_start",stage,id,debug)
    af.printl("S\t"+fasta_n_first+"\t"+str(range['start']),stage,id,debug)
    af.print_head_and_nlines(temp_mpileup, "mpileup file in reduce function", 10, stage, id, debug)
    af.file_and_folder_size("/tmp","before mpileup_merge_reducev3.sh", stage, id, debug)
    af.copy_to_s3(stage, id, debug, storage, BUCKET_NAME, temp_mpileup, temp_to_s3, s3_temp_folder)
    af.printl("Starting Merge Script",stage,id,debug)
    sinple_out = sp.check_output(['bash', '/function/bin/mpileup_merge_reducev3_nosinple.sh', temp_mpileup, '/function/bin/', buffer_size])
    sinple_out = sinple_out.decode('UTF-8')
    af.printl("Finished Merge Script",stage,id,debug)
    sinple_name=temp_mpileup+'_merged.mpileup'

    # write output to /tmp
    with open(sinple_name, 'w') as f:
        f.write(sinple_out)
    tr3 = af.execution_time("merge mpileups","time_end",stage,id,debug)
    af.printl("s\t"+fasta_n_first+"\t"+str(range['start']),stage,id,debug)
    af.printl(f' merge mpileups: execution_time:  {tr3 - tr2}: s',stage,id,debug)
    
    af.file_and_folder_size("/tmp","after mpileup_merge_reducev3.sh", stage, id, debug)
    # copy output to s3
    af.copy_to_s3(stage, id, debug, storage, BUCKET_NAME, sinple_name, temp_to_s3, s3_temp_folder)

    os.chdir(wd)
    #Upload part
    s3 = boto3.client('s3')
    part = s3.upload_part(
        Body = sinple_out,
        Bucket = BUCKET_NAME,
        Key = mpu_key,
        UploadId = mpu_id,
        PartNumber = n_part
    )
    tr4 = time.time()
    af.printl(f'reduce 2 - lambda: execution_time_total: {tr4 - tr0}: s',stage,id,debug)
    if (live_log==True):
        red1_log.stop()
    return {"PartNumber" : n_part, "ETag" : part["ETag"], "mpu_id": mpu_id}



if __name__ == "__main__":

    # Preliminary steps:
    # 1. upload the fastq.gz file(s) into BUCKET_NAME/fastqgz
    # 2. Upload the fasta file into BUCKET_NAME/fasta
    if (live_log==True):
        local_log = util.RemoteLoggingFeed(lambda_log_stream)  #
        local_log.start()
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
    ###################################################################
    # 0. Start the VM hosting the Redis Server
    ec2.action_on_instance(region_name=REGION_NAME,instance_id=EC2_ID,action='start')
    time.sleep(25)
    ec2.flush_redis_database(path_to_secret_key=PEM,ec2_IP_address=EC2_IP,redis_password='lucio-redis-server')
    time.sleep(5)


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
    fasta_chunks_prefix = re.sub("\.fasta|\.fa|\.fas", "_" + str(fasta_chunk_size) + "split_", fasta_file)
    fasta_list = lithopsgenetics.prepare_fasta(cloud_adr, runtime_id, FASTA_BUCKET, fasta_folder, fasta_file, split_fasta_folder, fasta_chunk_size, fasta_chunks_prefix, fasta_char_overlap)
    #fasta_metadata = storage.head_object('bucket', 'key')
    #total_fasta_size = int(fasta_metadata ['content-length'])
    print("Fasta file size: ")
    t3 = af.execution_time("fasta list","time_end",stage,id,debug)
    
    print(f'PP:0: fasta list: execution_time: {t3 - t2}: s')

    ###################################################################
    ###################################################################
    # 3. GENERATE ITERDATA
    print("\nGenerating iterdata")
    t4 = af.execution_time("iterdata","time_start",stage,id,debug)
    iterdata = lithopsgenetics.generate_alignment_iterdata(fastq_list, fasta_list, iterdata_n)
    
    
    # 4. SPLIT ITERDATA IF NUMBER HIGHER THAN CONCURRENT FUNCTIONS
    iterdata_n=len(iterdata)
    fastq_set_n=len(fasta_list) # number of fastq files aligned to each fasta chunk.
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
    ###################################################################
    # 4. PREPROCESSING SUMMARY
    

    #print(json.dumps(iterdata, indent=4))
    print("\nITERDATA LIST")
    #print(str(iterdata))
    iter_count=0
    for el in iterdata:
        iter_count+=1
        #print(str(iter_count)+":\t"+str(el))
    print("number of fasta chunks: " + str(fastq_set_n))
    print("number of fastq chunks: " + str(len(fastq_list)))
    print("fasta x fastq chunks: "+ str(len(fasta_list)*len(fastq_list)))
    print("number of iterdata elements: " + str(iterdata_n))
    PP_end = af.execution_time("preprocessing","time_start",stage,id,debug)
    print(f'PP:0: preprocessing: execution_time: {PP_end - PP_start}: s')

    if pre_processing_only==False:

        #ec2.start_redis_reducer(path_to_secret_key='~/bioss/cloudbutton/variant_caller/v3/aws-redis-damien.pem',ec2_IP_address=EC2_IP,iteration_number=iterdata_n)
        #time.sleep(5)

        ###################################################################
        ###################################################################
        # 5. MAP-REDUCE
        stage="MR"
        id=0
        start = af.execution_time("map-reduce","time_start",stage,id,debug)
        mapreduce = MapReduce(map_alignment, reduce_function, create_sinple_key, runtime_id, runtime_mem, runtime_mem_r, buffer_size, file_format, func_timeout_map, func_timeout_reduce, 'DEBUG', BUCKET_NAME, stage, id, debug, skip_map, EC2_IP,10000,lb_method,PEM, fastq_set_n, run_id, iterdata_n, redis_index_correction)
        map_time, creating_keys_time, reduce_time, multipart_upload_time = mapreduce(iterdata, iterdata_sets)
        end = af.execution_time("map-reduce","time_end",stage,id,debug)

        ###################################################################
        ###################################################################
        # 6. MAP-REDUCE SUMMARY AND REDIS LOG

        print("MAP-REDUCE SUMMARY")
    
        af.printl("map phase: execution_time_total_varcall: "+str(map_time)+": s",stage, id, debug)
        af.printl("reduce 1 phase: execution_time_total_varcall: "+str(creating_keys_time)+": s",stage, id, debug)
        af.printl("reduce 2 phase: execution_time_total_varcall: "+str(reduce_time)+": s",stage, id, debug)
        af.printl("reduce 3 phase: execution_time_total_varcall: "+str(multipart_upload_time)+": s",stage, id, debug)
        af.printl("total time: execution_time_total_varcall: "+str(end - start)+": s",stage, id, debug)
    
        # print redis log
        # r = Storage(backend='redis')
        # # af.printl("REDIS LOG")
        # # af.printl(r.get_object('',"redis_log.txt"))
        # # af.printl("END OF REDIS LOG")
        # redis_log = r.get_object('',"redis_log.txt")
        # af.print_formatted("REDIS LOG", redis_log)

        ###################################################################
        ###################################################################

    # CLEAN AND STOP THE VM HOSTING THE REDIS SERVER
    ec2.flush_redis_database(path_to_secret_key=PEM, ec2_IP_address=EC2_IP,redis_password='lucio-redis-server')
    time.sleep(5)
    ec2.action_on_instance(region_name=REGION_NAME,instance_id=EC2_ID,action='stop')
    if (live_log==True):
        local_log.stop()
