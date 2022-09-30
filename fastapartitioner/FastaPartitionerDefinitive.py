
import re
import lithops
import json
import gc
import time






class PartitionFasta():

        


    def __init__(self, storage, my_bucket_name,obj):
        self.storage = storage
        self.my_bucket_name = my_bucket_name
        self.obj = obj
        

    #Generate metadata from fasta file
    def generate_chunks(self,obj,total_obj_size,obj_size,overlap):
        data = obj.data_stream.sb.read().decode('utf-8')
        FASTA_HEADER = r">.+\n"
        content = ""
        found = False
        titles = list(re.finditer(FASTA_HEADER,data))
        i = 0
        start = 0
        last_el = 0
        prev = ""
        for m in titles:
            found = True
            start =  obj.data_byte_range[0] + m.start()
            if i == 0:
                if obj.part >= 1:
                    if start > (obj.part)*obj.chunk_size: 
                        content =  content +  str((obj.part)*obj.chunk_size) +',' + str(start+overlap)+ '\n'
            if (obj.part == int(total_obj_size)) and  (int(obj_size) - int(start + len(m.group()))) < total_obj_size:
                content =  content + str(start) +',' + str(obj_size)+ 'x\n'
            else:
                if str(prev) != str(start) and prev != "":
                    content =  content + prev +',' + str(start)+ '\n'
                    content =  content + str(start) +',' + str(start + len(m.group()))+ 'x\n'
                else:
                    content =  content + str(start) +',' + str(start + len(m.group()))+ 'x\n'
                prev = str(start + len(m.group()))

                
                
            last_el = start + len(m.group())
            i=+1
        
        if last_el < (obj.part+1)*obj.chunk_size and found:
            if obj.part == int(total_obj_size):
                content = content + str(last_el) +','+ str(int(obj_size))+ '\n'
            else:
                content = content + str(last_el) +','+ str(str((obj.part+1)*obj.chunk_size + overlap))+ '\n'

        if not found :
            if obj.part == int(total_obj_size):
                content = content + str(obj.data_byte_range[0]) +',' + str(obj_size)+ '\n'
            else:
                content = content + str(obj.data_byte_range[0]) +','+ str(((obj.part+1)*obj.chunk_size) + overlap)+ '\n'
        

        self.storage.put_object(self.my_bucket_name, f'cache/obj{obj.part}.data',content)
        
        return f'cache/obj{obj.part}.data'
       


def print_chunk_info(obj):
    # obj is a CloudObject
    print("===================")
    print(obj)
    print(obj.key)
    print('part:'+str(obj.part))
    print(obj.data_byte_range)  
    print(obj.chunk_size)  # in bytes
    print("===================")

def run_worker_metadata(obj,storage,my_bucket_name,total_obj_size,obj_size):
    print_chunk_info(obj)
    partitioner = PartitionFasta(storage,my_bucket_name,obj)
    dades = partitioner.generate_chunks(obj,total_obj_size,obj_size)
    return dades




    
#s3://sra-pub-src-2/SRR8774337/cleaned.fasta

if __name__ == "__main__":
    fexec = lithops.FunctionExecutor(
        log_level='DEBUG',runtime='testcontainer',runtime_memory=4096, max_workers=1000)
    
    obj = 's3://ayman-lithops-meta-cloudbutton-hutton/fasta/hg19.fa'  # Change-me
    my_bucket_name = 'ayman-lithops-meta-cloudbutton-hutton' # Change-me
    
    location = obj.split('//')[1].split('/')
    for_head = location[1:]
    for_head = '/'.join(for_head)
    data_bucket_name = location[0]  
    ovlp = 300
    storage = lithops.Storage()
    worker_chunk_size = 500000000
    fasta = storage.head_object(data_bucket_name, for_head) 
    total_fasta_size = int(fasta['content-length'])/worker_chunk_size
    seq_name = location[-1].split('.')[0]
  
    
    print('===================================================================================')
    print('# metadata chunks: '+ str(int(total_fasta_size)))
    print('bucket to access data: ' + str(data_bucket_name))
    print('reference genome name: ' + str(seq_name))
    print('fasta size: ' + str(fasta['content-length']) + ' bytes')
    print('===================================================================================')

    fexec.map(run_worker_metadata, {'my_bucket_name':my_bucket_name,'obj':obj,'total_obj_size':total_fasta_size,'obj_size':fasta['content-length']},
            obj_chunk_size=worker_chunk_size)
    fexec.wait()
    results = fexec.get_result()

    for el in results:
        for k in el.splitlines():
            if 'x' in k:
                print(k)

    fexec.plot()
    fexec.clean()
