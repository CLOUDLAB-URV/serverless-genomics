#
# (C) Copyright Cloudlab URV 2021
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

import pathlib
import os
import shutil
import tempfile
from lithops import Storage
import subprocess as sp
import lithopsgenetics.auxiliaryfunctions as af
import re
import Metadata

CWD = os.getcwd()
CURRENT_PATH = str(pathlib.Path(__file__).parent.resolve())
LOCAL_TMP = os.path.realpath(tempfile.gettempdir())


def preprocess_gzfile(cloud_adr, bucket_name, fastq_folder, index_folder, file_name):
    """
    The function takes the gzip file, creates the necessary index files
    for partitioning and stores them in the bucket.
    """
    print("PREPROCESSING FASTQ.GZ FILE")
    print("generating .gz index with gztool and uploading index to cloud storage")
    print("start of *preprocess_gzfile* function")
    os.chdir(LOCAL_TMP)
    storage = Storage()
    # 0 DOWNLOAD FASTQ.GZ FILE TO LOCAL TMP
    local_filename = os.path.join(LOCAL_TMP, file_name)
    if not os.path.isfile(local_filename):
        print(f'Downloading {cloud_adr}://{bucket_name}/{fastq_folder}/{file_name} to {LOCAL_TMP}')
        obj_stream = storage.get_object(bucket_name, fastq_folder + file_name, stream=True)
        with open(local_filename, 'wb') as fl:
            shutil.copyfileobj(obj_stream, fl)

    # 1 GENERATING THE INDEX AND INFORMATION FILES AND UPLOADING TO THE BUCKET
    sp.run(CURRENT_PATH+'/generateIndexInfo.sh '+local_filename, shell=True, check=True, universal_newlines=True)
    output = sp.getoutput(CURRENT_PATH+'/generateIndexInfo.sh '+local_filename)
    output = output.split()
    total_lines = str(af.only_numerics(output[-3]))

    # 2. UPLOAD FILES TO THE BUCKET
    remote_filename = index_folder + file_name

    print(f'Uploading {local_filename}i to {cloud_adr}://{bucket_name}/{index_folder}')
    with open(f'{local_filename}i', 'rb') as fl:
        storage.put_object(bucket_name, f'{remote_filename}i', fl)

    print(f'Uploading {local_filename}i.info to {cloud_adr}://{bucket_name}/{index_folder}')
    with open(f'{local_filename}i.info', 'rb') as fl:
        storage.put_object(bucket_name, f'{remote_filename}i.info', fl)

    print(f'Uploading {local_filename}i_tab.info to {cloud_adr}://{bucket_name}/{index_folder}')
    with open(f'{local_filename}i_tab.info', 'rb') as fl:
        storage.put_object(bucket_name, f'{remote_filename}i_tab.info', fl)

    os.remove(f'{local_filename}i')
    os.remove(f'{local_filename}i.info')
    os.remove(f'{local_filename}i_tab.info')

    os.chdir(CWD)
    print("end of *preprocess_gzfile* function")
    return total_lines


def download_index_files(bucket_name, index_folder, file_name, storage=None):
    """
    Downloads index files of a given fastq file if not present in the local machine
    """
    os.chdir(LOCAL_TMP)

    storage = Storage() if not storage else storage
    print(f'\n--> Downloading {file_name} index files')

    remote_filename = index_folder+file_name

    if not os.path.isfile(f'{file_name}i'):
        remote_filename_i = f'{remote_filename}i'
        obj_stream = storage.get_object(bucket_name, remote_filename_i, stream=True)
        with open(f'{file_name}i', 'wb') as fl:
            shutil.copyfileobj(obj_stream, fl)

    if not os.path.isfile(f'{file_name}i.info'):
        remote_filename_i = f'{remote_filename}i.info'
        obj_stream = storage.get_object(bucket_name, remote_filename_i, stream=True)
        with open(f'{file_name}i.info', 'wb') as fl:
            shutil.copyfileobj(obj_stream, fl)

    if not os.path.isfile(f'{file_name}i_tab.info'):
        remote_filename_i = f'{remote_filename}i_tab.info'
        obj_stream = storage.get_object(bucket_name, remote_filename_i, stream=True)
        with open(f'{file_name}i_tab.info', 'wb') as fl:
            shutil.copyfileobj(obj_stream, fl)

    os.chdir(CWD)


def split_fastafile(cloud_adr, runtime, bucket_name, fasta_folder, fasta_file, split_fasta_folder, fasta_chunk_size, chunk_overlap):
    """
    Splits fasta file into chunks depending on the desired chunk size
    """
    storage = Storage()

    print("start of *split_fastafile* function")
    print("--> Splitting fasta file")
    dir_list = os.listdir(".")


    # prints all files
    args = "python FastaPartitioner.py --data_location "+"s3"+"://"+bucket_name+"/"+fasta_folder+fasta_file+" --chunk_size "+str(fasta_chunk_size)+" --overlap "+str(chunk_overlap)+" --mybucket "+bucket_name+" --fasta_folder "+split_fasta_folder+" --runtime "+runtime
    args = args.split(" ")

    total_objects = sp.check_output(args).decode('utf-8')

    i = 0

    storage = Storage()
    data = ""
    chr_title = ""
    ini = 1
    while i <= int(total_objects):

        ret = storage.get_object(bucket_name, f'cache/obj{i}.data', stream = True).read().decode('utf-8')
        for line in ret.splitlines():
            if line == '\n' or line == '':
                pass
            elif '-' in line and ',' in line:
                d = line.split(',')
                range = d[1].split('-')
                abs_pos = int(range[1]) -  int(range[0])+1
                range = '\t'.join(range)
                data =   data + d[0].replace('>','')+ '\t' +chr_title[1].replace('>','') + '\t' + str(ini) + '\t' +  str(abs_pos)  + '\n'
                ini = abs_pos
            elif '-' in line:
                range = line.split('-')
                abs_pos = int(range[1]) -  int(range[0])+1
                data =  data + chr_title[0].replace('>','') + '\t' + chr_title[1].replace('>','') + '\t' + str(ini) + '\t' +  str(abs_pos) +'\n'
                ini = abs_pos
            else:
                chr_title = line.split(',')
                ini = 1
        i= i+1

    print(data)
    storage.put_object(bucket_name, f'cache/hashtable.data',data)
    print("end of *split_fasta* function")



def generate_fastq_chunk_list(bucket_name, index_folder, fastq_file_name, fastq_lines, seq_type, fastq_file_name2=None):
    """
    Creates the fastq chunk list (byte ranges)
    """
    print("start of *generate_fastq_chunk_list* function")
    print("CREATING FASTQ CHUNK BYTE RANGES LIST FROM FASTQ.GZ INFO FILES")
    storage = Storage()

    # 1 DOWNLOAD FASTQ INDEX FILES
    download_index_files(bucket_name, index_folder, fastq_file_name, storage)
    if seq_type == "paired-end":
        download_index_files(bucket_name, index_folder, fastq_file_name2, storage)

    os.chdir(LOCAL_TMP)
    total_lines = af.get_total_lines(fastq_file_name)
    total_lines = str(int(total_lines) + 1)
    print("\ntotal lines of fastq1 file: " + total_lines)
    total_lines2 = ""
    if seq_type == "paired-end":
        total_lines2 = af.get_total_lines(fastq_file_name2)
        total_lines2 = str(int(total_lines2) + 1)
        print("total lines of fastq2 file: " + total_lines2)
        if total_lines != total_lines2:
            raise Exception("paired end fastq line counts do not match: " + str(total_lines) + " vs " + str(total_lines2))

    # 2 GENERATE LINE INTERVAL LIST AND GET CHUNK BYTE RANGES
    print('\n--> Generating chunks (generateChunks.sh)')
    block_length = str(fastq_lines)
    sp.run('/home/agabriel/workspace/genomics_serverless/lithopsgenetics/generateChunks.sh '+fastq_file_name+' '+block_length+' '+total_lines, shell=True, check=True, universal_newlines=True)
    print('\n--> fastq 1 chunks information (af.read_chunks_info)')
    chunks, chunk_counter = af.read_chunks_info(fastq_file_name)
    print("fastq 1 chunks function output - chunk no.: " +  str(len(chunks)))
    print("fastq 1 chunks function output: " +  str(chunks))
    # processing second fastq file if present
    chunks2 = ""
    chunk_counter = ""
    if seq_type == "paired-end":
        # fastq file 2
        sp.run(CURRENT_PATH+'/generateChunks.sh '+fastq_file_name2+' '+block_length+' '+total_lines2, shell=True, check=True, universal_newlines=True)
        print('\n--> fastq 2 chunks information (af.read_chunks_info)')
        chunks2, chunk_counter2 = af.read_chunks_info(fastq_file_name2)
        print("fastq 2 chunks function output - chunk no.: " +  str(len(chunks2)))
        print("fastq 2 chunks function output: " +  str(chunks2))

    list_fastq = []
    for chunk in chunks:
        list_fastq.append((fastq_file_name, chunk))
    af.printl("fastq 1 list for iterdata: " + str(list_fastq))
    if seq_type == "paired-end":
        # fastq file 2 chunk list
        list_fastq2 = []
        for chunk2 in chunks2:
            list_fastq2.append((fastq_file_name2, chunk2))
        af.printl("fastq 2 list for iterdata: " + str(list_fastq2))
        list_fastq = [(i , j) for i, j in zip(list_fastq, list_fastq2)]
        af.printl("fastq 1 + 2 list for iterdata: " + str(list_fastq))

    return list_fastq


def generate_alignment_iterdata(list_fastq, list_fasta, iterdata_n):
    """
    Creates the lithops iterdata from the fasta and fastq chunk lists
    """
    os.chdir(CWD)
    iterdata = []
    
    # Number of fastq chunks processed. If doing a partial execution, iterdata_n will need to be multiple of the number of fasta chunks
    # so the index correction is done properly.
    num_chunks = 0  
    
    # Generate iterdata
    for fastq_key in list_fastq:
        num_chunks += 1
        for fasta_key in list_fasta:
            iterdata.append({'fasta_chunk': fasta_key, 'fastq_chunk': fastq_key})

    # Limit the length of iterdata if iterdata_n is not null.
    if iterdata_n is not None: 
        iterdata = iterdata[0:int(iterdata_n)]
        if(len(iterdata)%len(list_fasta)!=0):
            raise Exception("Hola")
        else: 
            num_chunks = len(iterdata)//len(list_fasta)

    return iterdata, num_chunks

def prepare_fastq(cloud_adr, BUCKET_NAME, idx_folder, fastq_folder, fastq_chunk_size, seq_type, fastq_file, seq_name,datasource,num_spots, fastq_file2=None):
    print("Creating fastq index files if they are not present")
    storage = Storage()
    index_status = ""

    try:
        if datasource == "s3":
            storage.head_object(BUCKET_NAME, idx_folder + fastq_file + "i")
    except:
        index_status = "fastq index not found"
    if datasource == "s3":
        if seq_type == "paired-end":
            try:
                storage.head_object(BUCKET_NAME, idx_folder + fastq_file2 + "i")
            except:
                index_status = "fastq index not found"
    elif datasource == "SRA":
        index_status = "fastq index not found"


    if index_status == "fastq index not found":
        if datasource == "s3":
            preprocess_gzfile(cloud_adr, BUCKET_NAME, fastq_folder, idx_folder, fastq_file)
        elif datasource == "SRA":
            chunks = preprocess_fastqsra(int(num_spots), int(fastq_chunk_size))
        # process second fastq file in the case of paired-end data
        if seq_type == "paired-end":
            if datasource == None:
                preprocess_gzfile(cloud_adr, BUCKET_NAME, fastq_folder, idx_folder, fastq_file2)
            elif datasource == "SRA":
                print("SRA as datasource; no need to preprocess the second fastq pair")
    else:
        index_status = "fastq index found"



    if datasource == "s3":
        print("fastq index status: " + index_status)
        fastq_list = generate_fastq_chunk_list(BUCKET_NAME, idx_folder, fastq_file, fastq_chunk_size, seq_type, fastq_file2)
        #print("fastq list for iterdata: " + str(fastq_list))
    elif datasource == "SRA":
        print("fastq index status: " + index_status)
        fastq_list = generate_fastq_chunk_list_fastq_sra(chunks,seq_name)
        #print("fastq list for iterdata: " + str(fastq_list))

    return fastq_list


def prepare_fasta(cloud_adr, runtime, BUCKET_NAME, fasta_folder, fasta_file, split_fasta_folder, fasta_chunk_size, fasta_chunks_prefix, fasta_line_overlap):
    storage = Storage()
    fasta_list = []
    try:
        print("Searching " + BUCKET_NAME + "/" + split_fasta_folder + fasta_chunks_prefix)
        print("Finding in storage: \n" + str(storage.list_keys(BUCKET_NAME, prefix=split_fasta_folder + fasta_chunks_prefix)))
        fasta_list = storage.list_keys(BUCKET_NAME, prefix=split_fasta_folder)
        #print("fasta list: " + str(fasta_list))
    except:
        print("split fasta folder empty / not found")


    if not fasta_list or fasta_list == []:
        print("splitting fasta file and storing chunks")
        split_fastafile(cloud_adr, runtime, BUCKET_NAME, fasta_folder, fasta_file, split_fasta_folder, fasta_chunk_size, fasta_line_overlap)
        try:
            print("Finding newly created fasta chunks")
            print("Searching " + BUCKET_NAME + "/" + split_fasta_folder + fasta_chunks_prefix)
            print("Finding in storage: \n" + str(storage.list_keys(BUCKET_NAME, prefix=split_fasta_folder + fasta_chunks_prefix)))
            fasta_list = storage.list_keys(BUCKET_NAME, prefix=split_fasta_folder + fasta_chunks_prefix)
            #print("fasta list: " + str(fasta_list))
        except:
            print("error generating fasta chunks")
    return fasta_list

#For fastq-dump and gztool
def fastq_to_mapfun(fastq_n, fastq_file_key, fastq_chunk_data, BUCKET_NAME, fastq_folder, idx_folder, datasource,stage,id,debug):
    '''
    Function executed within the map function to retrieve the relevant fastq chunk from object storage
    '''
    
    storage = Storage()
    if datasource == "s3":
        af.printl(fastq_n + " file key: " + str(fastq_file_key),stage,id,debug)
        af.printl(fastq_n + " chunk data: " + str(fastq_chunk_data),stage,id,debug)
        byte_range = f"{int(fastq_chunk_data['start_byte'])-1}-{int(fastq_chunk_data['end_byte'])}"
        temp_fastq_gz = af.copy_to_runtime(storage, BUCKET_NAME, fastq_folder, fastq_file_key, byte_range, id=id, debug=debug)

        # getting index and decompressing fastq chunk
        temp_fastq = temp_fastq_gz.replace('.fastq.gz', f'_chunk{fastq_chunk_data["number"]}.fastq')
        temp_fastq_i = af.copy_to_runtime(storage, BUCKET_NAME, idx_folder, f'{fastq_file_key}i', id=id, debug=debug, stage=stage)
        block_length = str(int(fastq_chunk_data['end_line']) - int(fastq_chunk_data['start_line']) + 1)
        cmd = f'gztool -I {temp_fastq_i} -n {fastq_chunk_data["start_byte"]} -L {fastq_chunk_data["start_line"]} {temp_fastq_gz} | head -{block_length} > {temp_fastq}'
        sp.run(cmd, shell=True, check=True, universal_newlines=True)
    elif datasource == "SRA":
        size = len(fastq_file_key)
        seq_name = fastq_file_key

        sp.call(['chmod','+x','fastq-dump'])
        # To supress the vdb-config warning
        sp.run(['vdb-config', '-i'])

        # Report cloud identity so it can take data from s3 needed to be executed only once per vm
        output = str(sp.run(['vdb-config', '--report-cloud-identity', 'yes'], capture_output=True).stdout)
        
        wd = os.getcwd()
        
        os.chdir(f"/tmp")
        temp_fastq = f'/tmp/'+seq_name+f'_chunk{fastq_chunk_data["number"]}.fastq'
        data_output = sp.run(['fastq-dump', str(seq_name), '-X', str(int(fastq_chunk_data["start_line"])) , '-N', str(int(fastq_chunk_data["end_line"])),'-O',f'/tmp'],
                             capture_output=True)
        #data_output = sp.run(['fastq-dump', str(seq_name), '-X', str(int(fastq_chunk_data["start_line"])) , '-N', str(int(fastq_chunk_data["end_line"])),'-Z'],
        #                     capture_output=True)
        af.printl("data_output contents: " + str(data_output),stage,id,debug)                     
        os.rename(f'/tmp/'+seq_name+'.fastq', temp_fastq)
        #af.printl('\n'.join(sorted(os.listdir("/tmp"))))
    
    return temp_fastq


#For Fastq-dump
def preprocess_fastqsra( num_spots, chunk_size):
    print("start of *preprocessfastq_file* function")

    end = 0
    ini = 0

    mod =  num_spots % chunk_size
    data = ""
    while end < num_spots:
        if end == num_spots-(chunk_size+mod):
            ini = end
            end = num_spots
            data = data + str(ini)+" "+str(end) + '\n'
        else:
            ini = end
            end = end + chunk_size
            data = data + str(ini)+" "+str(end) + '\n'

    print("end of *preprocessfastq_file* function")
    return data

#For fastq-dump
def generate_fastq_chunk_list_fastq_sra(data,seq_name):
    """
    Creates the fastq chunk list from a file
    """
    chunks, chunk_counter = af.read_chunks_info_random(data)
    #print("fastq chunks function output - chunk no.: " + str(chunk_counter))
    #print("fastq chunks function output: " + str(chunks))

    list_fastq = []
    for chunk in chunks:
        #print((f'{seq_name}' , chunk))
        list_fastq.append((seq_name , chunk))

    return list_fastq
