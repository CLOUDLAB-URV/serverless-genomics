import subprocess as sp
import os
from typing import List

def prepare_fastq(fastq_chunk_size: int, seq_name: str, num_spots: int) -> List[str]:
    """
    Prepare fastq chunks for processing. The fastq files can be obtained from SRA.
    """
    chunks = preprocess_fastqsra(int(num_spots), int(fastq_chunk_size))
    fastq_list = generate_fastq_chunk_list_fastq_sra(chunks,seq_name)

    return fastq_list


def preprocess_fastqsra(num_spots: int, chunk_size: int) -> str:
    """
    Get fastq data.
    """
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
    return data


def generate_fastq_chunk_list_fastq_sra(data: str, seq_name: str) -> List[str]:
    """
    Creates the fastq chunk list
    """
    chunks = read_chunks_info_random(data)

    list_fastq = []
    for chunk in chunks:
        list_fastq.append((seq_name , chunk))

    return list_fastq


def read_chunks_info_random(data: str) -> List[str]:
    chunk_list = []
    counter = 0
    
    for line in data.splitlines():
        info = line.split(' ')
        chunk_list.append({'number': (counter+1), 'start_line': str(info[0]),
                            'end_line': str(info[1])})
        counter += 1
    return chunk_list


def fastq_to_mapfun(fastq_file_key: str, fastq_chunk_data: str) -> str:
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