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

    return last_line.split(' ')[1]


# Function to read information of the chunks from a file
def read_chunks_info(file):
    chunk_list = []
    counter = 0
    with open(file+'.chunks.info', 'r') as f:
        for line in f:
            #print('Chunk info:'+str(counter))
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


def copy_to_runtime(storage, bucket, folder, file_name, byte_range=None):
    print(f'Copying {file_name} to local disk')
    extra_get_args = {'Range': f'bytes={byte_range}'} if byte_range else {}
    obj_stream = storage.get_object(bucket=bucket, key=folder+file_name, stream=True, extra_get_args=extra_get_args)
    temp_file = "/tmp/" + file_name
    with open(temp_file, 'wb') as file:
        shutil.copyfileobj(obj_stream, file)
        print(f'Finished copying {file_name} to local disk')
    return temp_file

def check_temp_file(temp_file, no_of_lines):
    print("\nPrinting " + temp_file )
    file = open(temp_file, 'r')
    Lines = file.readlines()
    count = 0
    # Strips the newline character
    for line in Lines:
        count += 1
        if count < no_of_lines:
            print(line.strip())
            #print("Line{}: {}".format(count, line.strip()))
    print("Finished printing " + temp_file + "\n")

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

def print_head_and_size(file_name, file_desc, number_of_lines):
    file = open(file_name)
    #number_of_lines = 5
    print("printing " + file_desc)
    print("printing " + file_name)
    size = os.path.getsize(file_name)
    print("file size " + str(size))
    for i in range(number_of_lines):
        line = file.readline()
        line = line.rstrip("\n")
        print(str(i) + "\t" + line)

def print_sp_output(header, file):
    print("\n"+ header) 
    file_list = str(file).split("\\n")
    for line in file_list:
        print(line)
    print("############## " + header + " END")

#Temporal function to test.
# Function to read information of the chunks from a file, without byteranges.
def read_chunks_info_random(file):
    chunk_list = []
    counter = 0
    with open(file, 'r') as f:
        for line in f:
            #print('Chunk info:'+str(counter))
            info = line.split()
            chunk_list.append({'number': (counter+1), 'start_line': str(info[0]),
                               'end_line': str(info[1])})
            counter += 1
    return chunk_list, counter

def delete_temp_files(dir_name):
    tmp_dir = os.listdir(dir_name)
    for item in tmp_dir:
        if item.endswith((".gz",".gzi",".fastq",".map",".txt",".fasta",".fai",".gem",".info")):
            os.remove(os.path.join(dir_name, item))