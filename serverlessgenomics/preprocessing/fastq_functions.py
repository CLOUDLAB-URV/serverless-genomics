import subprocess as sp
import os
from typing import List

def prepare_fastq(chunk_size: int, seq_name: str, num_spots: int) -> List[str]:
    """
    Prepare fastq chunks for processing. The fastq files can be obtained from SRA.
    """
    list_fastq = []
    ini = end = counter = 0
    mod = num_spots % chunk_size

    while end < num_spots:
        ini = end
        end = num_spots if end == num_spots - (chunk_size + mod) else end + chunk_size
        list_fastq.append(( seq_name, {'number': (counter+1), 'start_line': str(ini), 'end_line': str(end)}))
        counter += 1
    return list_fastq
  

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