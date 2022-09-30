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
        
        
    
    def get_object_size(self,file_name):
        obj_metadata = self.storage.head_object(self.bucket_name, file_name)
        return obj_metadata['content-length']

    def utf8len(self,s):
        return len(s.encode('utf-8'))

    #Generate metadata from fasta file
    def generate_metadata(self,element,range):
        data = self.storage.get_object(self.bucket_name, element, extra_get_args={ 'Range': 'bytes='+ str(range['start']) + '-' + str(range['end']) },stream=True)
        print(len(data))
        print(data)
        FASTA_HEADER = r"^>.*"
        byte_range_end = range['start']
        print(range)
        
        file = ""
        bytes_part = range['start']
        for line in data.splitlines(keepends=True):
            byte_range_end = byte_range_end + self.utf8len(line) 
            found = re.search(FASTA_HEADER,line)
            if found:
                print("if")
                file = file + '\"bytes_part\"' + ':' + str(bytes_part) + ',\n'
                file =  file + '\"title\":'+' {\"start\":'+ str(bytes_part) +','+ ' \"end\":' + str(byte_range_end)+ '}\n'
                bytes_part = byte_range_end
            else:
                print('else')
                bytes_part= bytes_part+self.utf8len(line)
            
            file = file + '\"bytes_part\"' + ':' + str(bytes_part) + ',\n'
        part = range['part']
        self.storage.put_object(self.bucket_name, 'fasta/'+str(part)+'.data',file)

    def reduce_metadata(self):
        pass
    
    def read_chunks(self):
        #Might generate memory errors, since it's reading the whole file at once, if the runtime memory is not enough, then it should be split into chunks
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
        self.storage.put_object(self.bucket_name, 'fasta/metadata.data',file.read())

        file.close()
        
    def read_metadata(self):
        import json
        f = open("metadata.json")
        data = json.load(f)
        return data
    
    def generate_chunks(self,metadata,chunk_size):
        total_chr = len(metadata)
        actual_chr = 0
        consumed_chromosomes = 0
        chunks = []
        row = []
        chunk_start = metadata[actual_chr]['title']['end']
        chunk_end = chunk_size
        chromosome_subindex = 0
        while consumed_chromosomes < total_chr:

            if metadata[actual_chr]['sequence']['end'] < chunk_end: 
                if (row == [] and chunks == []):
                    dict = {
                        'title': 
                            {
                                'start':metadata[actual_chr]['title']['start'], 
                                'end':metadata[actual_chr]['title']['end'] 
                            },
                        'sequence': 
                            {
                                'start':metadata[actual_chr]['title']['end'], 
                                'end':metadata[actual_chr]['sequence']['end']
                            },
                            'chromosome_index': actual_chr,
                            'chromosome_subindex': chromosome_subindex
                    }
                    row.append(dict)
                    actual_chr = actual_chr + 1
                    chromosome_subindex = chromosome_subindex + 1
                    consumed_chromosomes = consumed_chromosomes + 1
                    if actual_chr < total_chr:
                        chunk_start = metadata[actual_chr]['title']['end']
                else:
                    if row != []:
                        dict = {
                        'title': 
                            {
                                'start':metadata[actual_chr]['title']['start'], 
                                'end':metadata[actual_chr]['title']['end'] 
                            },
                        'sequence': 
                            {
                                'start':metadata[actual_chr]['title']['end'], 
                                'end':metadata[actual_chr]['sequence']['end']
                            },
                            'chromosome_index': actual_chr,
                            'chromosome_subindex': chromosome_subindex
                        }
                        chromosome_subindex = chromosome_subindex + 1
                        row.append(dict)
                    else:
                        dict = {
                        'title': 
                            {
                                'start':metadata[actual_chr]['title']['start'], 
                                'end':metadata[actual_chr]['title']['end'] 
                            },
                        'sequence': 
                            {
                                'start':chunks[-1][-1]['sequence']['end'], 
                                'end':metadata[actual_chr]['sequence']['end']
                            },
                            'chromosome_index': actual_chr,
                            'chromosome_subindex': chromosome_subindex
                        }
                    chromosome_subindex = chromosome_subindex + 1
                    row.append(dict)
                    actual_chr = actual_chr + 1
                    consumed_chromosomes = consumed_chromosomes + 1
                    chromosome_subindex = 0
                    if actual_chr < total_chr:
                        chunk_start = metadata[actual_chr]['title']['end']
            else:
                dict = {
                        'title': 
                            {
                                'start':metadata[actual_chr]['title']['start'], 
                                'end':metadata[actual_chr]['title']['end'] 
                            },
                        'sequence': 
                            {
                                'start':chunk_start, 
                                'end':chunk_end
                            },
                            'chromosome_index': actual_chr,
                            'chromosome_subindex': chromosome_subindex

                        }
                chromosome_subindex = chromosome_subindex + 1
                chunk_start = chunk_end        
                chunk_end = chunk_end + chunk_size
                row.append(dict)

                chunks.append(row)
                
                row = []  
        
        chunks.append(row)
        
        i=0
        for chunk in chunks:
            
            print(chunk)
            i=i+1
        print('===================================================================================')
        print('===================================================================================')
        print('===================================================================================')
        print('===================================================================================')
        
        return chunks
        
    #Overlap only happens if a chromosome is separated into multiple files.
    def resize_chunks(self,chunks,overlap):
        i = 0
        chunk_counter = 0
        for chunk in chunks:
            chunk_size = len(chunk)-1
            for el in chunk:
                if chunk_size > 1 and i+1 < chunk_size:
                    if el['chromosome_index'] == chunk[i+1]['chromosome_index']:
                        el['sequence']['end'] = el['sequence']['end'] + overlap                        
                else:
                    j = 0
                    if(chunk_counter < len(chunks)-1):
                        for el2 in chunks[chunk_counter+1]:
                            if el2['chromosome_index'] == el['chromosome_index']:
                                el['sequence']['end'] = el2['sequence']['start'] + overlap
                                break
                            j = j + 1
                i=+1
            chunk_counter=+1
        i=0
        for chunk in chunks:
            print(chunk)
            i=i+1

   
                
                 
            
                    
#vim /tmp/lithops/logs/3eb5b2-0-A000.log      

def print_chunk_info(obj):
    # obj is a CloudObject
    print(obj)
    print(obj.key)
    print(obj.part)
    print(obj.data_byte_range)  
    print(obj.chunk_size)  # in bytes



def run_worker_metadata(obj,storage,bucket_name):
    print_chunk_info(obj)
    partitioner = PartitionFasta(storage,bucket_name,obj)
    partitioner.generate_metadata(obj)

def reduce_metadata(obj,storage,bucket_name):
    partitioner = PartitionFasta(storage,bucket_name,obj)
    partitioner.reduce_metadata()

def run_worker(obj, storage,bucket_name):
    partitioner = PartitionFasta(storage,bucket_name,obj)
    partitioner.read_chunks()
    metadata = partitioner.read_metadata()
    chunks = partitioner.generate_chunks(metadata,2000000)
    partitioner.resize_chunks(chunks,300)
    return 0

    
def get_object_size(obj,storage,bucket_name):
    partitioner = PartitionFasta(storage,bucket_name,obj)
    return partitioner.get_object_size(obj)
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
        log_level='DEBUG',runtime='testcontainer',storage='aws',backend='aws_lambda')
    
    args = {'obj':'/fasta/hg19.fa', 'bucket_name':'ayman-lithops-meta-cloudbutton-hutton'}
    
    iterdata = []
    
    file_size = int(get_object_size(args['obj'],fexec.storage,args['bucket_name']))
    
    start = 0
    end = 0
    i = 0
    while start < int(file_size):
        end = start + 2000000
        if end > file_size:
            end = file_size
        iterdata.append({'start':start,'end':end,'part':i})
        start = end
        i = i + 1
    

    fexec.map(run_worker,iterdata,extra_args=args)
    fexec.wait()
    res = fexec.get_result()
 
    fexec.clean()



