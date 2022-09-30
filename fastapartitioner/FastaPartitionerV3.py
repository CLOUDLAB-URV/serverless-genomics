import chunk
import re
import lithops
from lithops import Storage
import ast

class PartitionFasta():

    
    def __init__(self, storage, bucket_name,obj):
        self.storage = storage
        self.bucket_name = bucket_name
        self.obj = obj
        
    #Reads fasta chunks from the original file that is located in S3, based in 2 regexs and size per chunk and returns the chunks
    #def read_fasta_chunks(self,obj,prev_header,index,chunk_index)
    def utf8len(self,s):
        return len(s.encode('utf-8'))

    def read_chunks(self):
        #Might generate memory errors, since it's reading the whole file at once
        data = self.storage.get_object(self.bucket_name, self.obj.path, stream=True)   
        byte_range_start = 0
        byte_range_end = 0
        FASTA_HEADER = r"^>.*" 
        file = open('metadata.json', 'w+')

        chr_num = 0
        locations = []
        
        for line in data:
            line = line.decode("utf-8")
            byte_range_end = byte_range_end + self.utf8len(line)
            found = re.search(FASTA_HEADER,line)
            if found:
                if chr_num == 0:
                    file.write('[\n')
                    file.write('\t' + '\t{\"title\":'+' {\"start\":'+ str(byte_range_start) +','+ ' \"end\":' + str(byte_range_end)+ '},\n')
                    byte_range_start = byte_range_end
                    chr_num = chr_num + 1
                elif chr_num > 0:
                    seq_end = byte_range_end-self.utf8len(line)
                    file.write('\t\t' +  '\"sequence\":'+' {\"start\":'+ str(byte_range_start) +','+ ' \"end\":' + str(seq_end)+ '}},\n')
                    file.write('\t\t' +  '{\"title\":'+' {\"start\":'+ str(seq_end) +','+ ' \"end\":' + str(byte_range_end)+ '},\n')
                    byte_range_start = byte_range_end
                    chr_num = chr_num + 1
        
        file.write('\t\t' +  '\"sequence\":'+' {\"start\":'+ str(byte_range_start) +','+ ' \"end\":' + str(byte_range_end)+ '}}\n')

        file.write(']')
        self.storage.put_object(self.bucket_name, f'fasta/metadata.data',file.read())

        file.close()
        
    def read_metadata(self):
        import json
        f = open("metadata.json")
        data = json.load(f)
        return data
    
    def generate_chunks(self,metadata,chunk_size):
        total_chr = len(metadata)
        actual_chunk = 0
        actual_chr = 0
        total_size = metadata[total_chr-1]['sequence']['end']
        consumed_chromosomes = 0
        chunks = []
        row = []
        
        accumulated_size = 0
        chunk_counter = 0
        part = False
        row.append(metadata[actual_chr])
        while actual_chr < total_chr:
            chromosome_size = metadata[actual_chr]['sequence']['end'] - metadata[actual_chr]['title']['start']
            chunk_size_count = chromosome_size
            if (chunk_size_count <= chunk_size):
                if (chunks == [] and row == []):
                    row.append(metadata[actual_chr])
                    actual_chr = actual_chr + 1
                else:
                    while (chunk_size_count <= chunk_size):
                        prev_chunk_size_count = chunk_size_count
                        chunk_size_count = chunk_size_count + metadata[actual_chr]['sequence']['end'] - metadata[actual_chr]['title']['start']
                        if chunk_size_count >= chunk_size:
                            dict = {
                                            'title': 
                                                {
                                                    'start':metadata[actual_chr]['title']['start'], 
                                                    'end':metadata[actual_chr]['title']['end'] 
                                                },
                                            'sequence': 
                                                {
                                                    'start':metadata[actual_chr]['sequence']['start'], 
                                                    'end':metadata[actual_chr]['sequence']['start'] + (prev_chunk_size_count - chunk_size_count)
                                                }
                                            }
                            part = True
                            row.append(dict)
                            chunk_size_count = chunk_size
                            chunk_counter =+ 1 
                        else:
                            dict = {
                                            'title': 
                                                {
                                                    'start':metadata[actual_chr]['title']['start'], 
                                                    'end':metadata[actual_chr]['title']['end'] 
                                                },
                                            'sequence': 
                                                {
                                                    'start':metadata[actual_chr]['sequence']['start'], 
                                                    'end':metadata[actual_chr]['sequence']['end']
                                                }
                                            }
                        
                            row.append(dict)
                        chunks.append(row)
        print(chunks)
                    
                    
               

def print_chunk_info(obj):
    # obj is a CloudObject
    print(str("obj part :") + str(obj.part))




def run_worker(obj, storage,bucket_name):
    partitioner = PartitionFasta(storage,bucket_name,obj)
    partitioner.read_chunks()
    metadata = partitioner.read_metadata()
    partitioner.generate_chunks(metadata,2000000)
    return 0

    
    """
    #Mapper phase 1: Read data using lithops.utils.WrappedStreamingBodyPartition and chunk it into FASTA chunks based on 2 regex patterns.
    
    

    #Mapper phase 2: Read previous metadata written by each chunk,retrieve chunks and write them into a file.
    #Prevent race condition between mappers, use of s3 waiters.
    FASTA_HEADER = r"^>.*"
    prev_mapper = obj.part - 1
    #First fasta chunk is always written, not caring about consistency.
    if obj.part == 0:
        wsb.write_fasta_chunk(res,obj,'>ch_',1,1)
    else:
    #Else, get the previous chunk header
        prev_chunk = wsb.get_chunk_waiting(prev_mapper).read().decode("utf-8")
        prev_header = re.findall(FASTA_HEADER, prev_chunk,re.MULTILINE)
   
    
        
    #Mapper phase 3: write the chunk into s3
        print("prev_header"+str(prev_header))
        wsb.write_fasta_chunk(res,obj,prev_header.pop(),1,1)

    return res
"""

if __name__ == "__main__":
    fexec = lithops.FunctionExecutor(
        log_level='DEBUG',runtime='testcontainer',storage='localhost',backend='localhost',max_workers=1)
    
    args = {'obj':'/mnt/c/Users/ayman/Desktop/fasta.fasta', 'bucket_name':'data'}
    
    
    fexec.call_async(run_worker, args)
    #fexec.map(run_worker, args,
              #obj_chunk_size=int(2))
    fexec.wait()
    res = fexec.get_result()
 
    fexec.clean()

