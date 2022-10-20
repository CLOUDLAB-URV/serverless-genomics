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
import shutil
import os
import re
import time
import subprocess as sp


def gzfile_info():
    file = input('Indicate the name of the gz file to chunk:')
    lines = input('Indicate the number of lines to retrieve for a chunk:')
    return file, lines


# Function to get input information
def gzfile_info_random():
    file = input('Indicate the name of the gz file to chunk:')
    start_line = input('Indicate the start line of the chunk:')
    end_line = input('Indicate the end line of the chunk:')
    return file, start_line, end_line


# Function to remove all non-numeric characters from a string
def only_numerics(string):
    string_type = type(string)
    return string_type().join(filter(string_type.isdigit, string))


def get_total_lines(file):
    """
    gets the total number of lines from the index_tab.info file
    """
    with open(file+'i_tab.info', 'r') as f:
        for last_line in f:
            pass

    return last_line.split(' ')[0].replace('(', '')


# Function to read information of the chunks from a file
def read_chunks_info(file):
    chunk_list = []
    counter = 0
    with open(file+'.chunks.info', 'r') as f:
        for line in f:
            printl('Chunk info:'+str(counter))
            info = line.split()
            chunk_list.append({'number': (counter+1), 'start_line': str(info[0]),
                               'end_line': str(info[1]), 'start_byte': str(info[2]),
                               'end_byte': str(info[3])})
            counter += 1
    return chunk_list, counter


# Function to read information of the chunks from a file
def read_chunk_info_random(file, start_line, end_line):
    with open(file+'.random_chunk_'+start_line+'_'+end_line+'.info', 'r') as f:
        for line in f:
            info = line.split()
            chunk = {'number': (1), 'start_line': str(info[0]), 'end_line': str(info[1]),
                     'start_byte': str(info[2]), 'end_byte': str(info[3])}
    return chunk

def execution_time(substage,start_end, stage, id, debug):
    timestamp=time.time()
    printl(substage+ ": "+start_end+":"+str(timestamp), stage, id, debug)
    return timestamp


def copy_to_runtime(storage, bucket, folder, file_name, stage, id, debug, byte_range=None, key=None):
    printl(f'Copying {file_name} to /tmp folder',stage,id,debug)
    obj_stream = storage.get_object(bucket=bucket, key=folder+file_name, stream=True, extra_get_args=byte_range)
    # Modified, temp_file = "/tmp/" + file_name
    temp_file = "/tmp/" + (f"{key}_" if key is not None else '') + file_name
    with open(temp_file, 'wb') as file:
        shutil.copyfileobj(obj_stream, file)
        printl(f'Finished copying {temp_file} to /tmp folder',stage,id,debug)
    return temp_file

def copy_to_s3(stage, id, debug, storage, BUCKET_NAME, file_name, temp_to_s3, folder_name=""):
    if temp_to_s3==True:
        printl("\ncopying "+file_name+" to "+BUCKET_NAME,stage,id,debug)

        destination_key= folder_name + os.path.basename(file_name)
        printl("destination key: "+destination_key,stage,id,debug)

        with open(file_name, 'rb') as file:
            data = file.read()
            storage.put_object(BUCKET_NAME, destination_key, body=data)
        return destination_key
    else:
        return 0

def check_temp_file(stage, id, debug, temp_file, no_of_lines):
    printl("\nPrinting " + temp_file,stage,id,debug)
    file = open(temp_file, 'r')
    Lines = file.readlines()
    count = 0
    # Strips the newline character
    for line in Lines:
        count += 1
        if count < no_of_lines:
            printl(line.strip(),stage,id,debug)
            #printl("Line{}: {}".format(count, line.strip()))
    printl("Finished printing " + temp_file + "\n",stage,id,debug)

def parse_optional_arg(arg_opt1, arg_opt2, text2_opt1=None, text2_opt2=None):
    out1 = ""
    out2 = ""
    if arg_opt1 is not None:
        out1 = arg_opt1
        out2 = text2_opt1
    else:
        out1 = arg_opt2
        out2 = text2_opt2
    if out2 is not None:
        return(out1, out2)
    else:
        return(out1)


#Temporal function to test.
# Function to read information of the chunks from a file, without byteranges.
def read_chunks_info_random(data):
    chunk_list = []
    counter = 0
    
    #print(data)
    for line in data.splitlines():
        #printl('Chunk info:'+str(counter))
        info = line.split(' ')
        #print(info)
        chunk_list.append({'number': (counter+1), 'start_line': str(info[0]),
                            'end_line': str(info[1])})
        counter += 1
    return chunk_list, counter
    


def getObject(key):
    storage = Storage()
    with open('/tmp/file.mpileup', 'a') as f:
        mpileup = storage.get_object(BUCKET_NAME, key).decode('UTF-8')
        f.write(mpileup)
        del mpileup  
		
def print_formatted(header, file,stage="",id="",debug=True):
    printl(header,stage,id,debug) 
    file_list = str(file).split("\\n")
    for line in file_list:
        line = re.sub(r'\\t', '\t', line)
        printl(line,stage,id,debug)
    printl("############## " + header + " END",stage,id,debug)

def print_formatted_redis(header, file,stage="",id="",debug=True):
    printl(header,stage,id,debug) 
    file_list = str(file).split("\\n")
    for line in file_list:
        line = re.sub(r'\\t', '\t', line)
        printl("REDIS\t"+line,stage,id,debug)
    printl("############## " + header + " END",stage,id,debug)


def print_indexer_output(header, file,length, stage, id, debug):
    printl(header,stage,id,debug) 
    file_list = str(file).split("\\n")
    for line in file_list:
        line = re.sub(r'\\t', '\t', line)
        if length=="short":
            if re.search(r"(%)", line)==None:
                printl(line,stage,id,debug)
        else:
            printl(line,stage,id,debug)
    printl("############## " + header + " END", stage, id, debug)


def print_stdout(header, file, stage="", id="", debug=True):
    printl("\n"+ header,stage,id,debug) 
    for line in file.stdout.rstrip():
        printl(line,stage,id,debug)
    #printl((file.stdout).rstrip())
    printl("############## " + header + " END",stage,id,debug)

def printl(string,stage="",id="",debug=True):
    '''
    add flush=True to print statement to enable real-time printing to log; printl = print to log
    '''
    if debug == True:
        if stage =="":
            print(string, flush=True)
            #print("")
        else:
            print(stage+"\t"+str(id)+"\t"+string, flush=True)
            #print("")

        

def file_and_folder_size(dir,desc, stage, id, debug):
    printl("",stage,id,debug)
    printl("##########",stage,id,debug)
    printl(dir.upper()+" FOLDER AND FILE SIZE - "+desc.upper(),stage,id,debug)
    cmd = f'ls -l '
    out=((sp.run(cmd+dir, shell=True, check=True, universal_newlines=True, capture_output= True)).stdout)
    out=out.splitlines()
    for line in out:
        printl(line,stage,id,debug)
    cmd = f'du -sh '
    printl("/tmp folder size, "+desc+": "+(sp.run(cmd+dir, shell=True, check=True, universal_newlines=True, capture_output= True)).stdout,stage,id,debug)
    

def print_head_and_nlines(file_name, file_desc, number_of_lines, stage, id, debug):
    file = open(file_name)
    line_count = len(file.readlines())
    cmd = f'head -n '+str(number_of_lines)+' '

    printl("##########",stage,id,debug)
    printl("Overview of " + file_desc,stage,id,debug)   
    printl(file_name+" total number of lines " + str(line_count),stage,id,debug)
    printl(file_name+" printing first "+str(number_of_lines)+" lines",stage,id,debug)
    out=(sp.run(cmd+file_name, shell=True, check=True, universal_newlines=True, capture_output= True)).stdout
    out=out.splitlines()
    for line in out:
        printl(line,stage,id,debug)

def clear_tmp(stage, id, debug):
    dir_name=f'/tmp'
    if os.path.exists(dir_name):
        tmp_dir = os.listdir(dir_name)
        for item in tmp_dir:
            if item.endswith((".gz",".gzi",".fastq",".map",".txt",".fasta",".gem",".info",".mpileup", ".csv")):
                os.remove(os.path.join(dir_name, item))
        printl("\n##########",stage,id,debug)
        printl("/tmp directory cleared",stage,id,debug)

