
import re
import lithops
from operator import itemgetter
import os
import json
import gc

class PartitionFasta():

    
    def __init__(self, storage, bucket_name,obj):
        self.storage = storage
        self.bucket_name = bucket_name
        self.obj = obj
        
        
    
    
    def utf8len(self,s):
        return len(s.encode('utf-8'))

    #Generate metadata from fasta file
    def generate_metadata(self,obj):
        data = obj.data_stream.sb.read().decode('utf-8')
        print('len_data:'+str(len(data)))
        print('byte_range:' + str(obj.data_byte_range[1]-obj.data_byte_range[0]))
        print('obj_size:' + str(obj.chunk_size))
        
        FASTA_HEADER = r"^>.*"
        byte_range_end = obj.data_byte_range[0]
        file = ""
        bytes_part = obj.data_byte_range[0]
        for line in data.splitlines(keepends=True):
            byte_range_end = byte_range_end + self.utf8len(line) 
            found = re.search(FASTA_HEADER,line)
            if found:
                file =  file + '{\"start\":'+ str(bytes_part) +','+ ' \"end\":' + str(byte_range_end)+ '}\n'
                bytes_part = byte_range_end
            else:
                bytes_part= bytes_part+self.utf8len(line)        
        if file != "":
            self.storage.put_object(self.bucket_name, f'fasta/obj{obj.part}.data',file)
        

    def reduce_metadata(self):
        metadata = set()
        list = []
        for filename in os.listdir('/tmp/lithops/data/fasta/'):
            f = open(f'/tmp/lithops/data/fasta/{filename}', 'r')
            for dict in f:
                metadata.add(dict)
        
        for el in metadata:
            list.append(json.loads(el))

        sorted_list = sorted(list, key=lambda d: d['start'])
        final_sorted_list = []
        del metadata
        del list
        gc.collect() 
                
        i=0
        size_sorted = len(sorted_list)

        while i < size_sorted:
            if i+1 < size_sorted:
                dict = {
                    "title": {
                        "start": sorted_list[i]['start'], 
                        "end": sorted_list[i]['end']
                        },
		            "sequence": {
                        "start":sorted_list[i]['end'], 
                        "end":sorted_list[i+1]['start']
                        }
                    }
                final_sorted_list.append(dict)
                
            else:
                dict = {
                    "title": {
                        "start": sorted_list[i]['start'], 
                        "end": sorted_list[i]['end']
                        },
                    "sequence": {
                        "start":sorted_list[i]['end'], 
                        "end":int(self.storage.head_object(self.bucket_name, self.obj.path)['content-length'])
                        }
                    }
                final_sorted_list.append(dict)
            i=i+1
        self.storage.put_object(self.bucket_name,'reduce_meta/metadata.json',json.dumps(final_sorted_list,indent=2))
            
        for el in final_sorted_list:
            print(el)
    
        
    def read_metadata(self):
        data = self.storage.get_object(self.bucket_name, f'reduce_meta/metadata.json', stream=True)
        
        return data
    
    def generate_chunks(self,metadata,chunk_size):
        metadata = json.load(metadata)
        print(metadata)
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
            print(f'chunk #{i}')
            print(chunk)
            print( 'size:'+str(chunk[-1]['sequence']['end']-chunk[0]['sequence']['start']))
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
            print(f'chunk #{i}')
            print(chunk)
            i=i+1

   
                
                 
            
                    
#vim /tmp/lithops/logs/3eb5b2-0-A000.log      

def print_chunk_info(obj):
    # obj is a CloudObject
    print(obj)
    print(obj.key)
    print(f'part:{obj.part}')
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
    metadata = partitioner.read_metadata()
    chunks = partitioner.generate_chunks(metadata,1*pow(2,25))
    partitioner.resize_chunks(chunks,300)
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
        log_level='DEBUG',runtime='testcontainer',max_workers=1000)
    
    args = {'obj':'s3://ayman-lithops-meta-cloudbutton-hutton/fasta/hg19.fa', 'bucket_name':'data'}
    
    
    fexec.map(run_worker_metadata, args,
              obj_chunk_size=1*pow(2,25))
    fexec.wait()

    fexec.call_async(reduce_metadata, args)
    fexec.wait()
    
    fexec.call_async(run_worker, args)
    #fexec.map(run_worker, args,
              #obj_chunk_size=int(2))
    fexec.wait()
    res = fexec.get_result()
    fexec.plot()
    fexec.clean()
