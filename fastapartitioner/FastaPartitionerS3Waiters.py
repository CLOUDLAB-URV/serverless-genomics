import re
from xxlimited import new
import lithops
from lithops import Storage


class WrappedStreamingBodyPartitionFasta():

    
    def __init__(self, wsb: lithops.utils.WrappedStreamingBodyPartition, storage, bucket_name):
        self.wsb = wsb
        self.storage = storage
        self.bucket_name = bucket_name

    #Reads fasta chunks from the original file that is located in S3, based in 2 regexs and size per chunk and returns the chunks
    #def read_fasta_chunks(self,obj,prev_header,index,chunk_index)
    #Modify since its buggy.
    def read(self):
        data = self.wsb.sb.read()      
        FIND_FASTA_HEADER = b">.+\n[acgtnACGTN\n]+|[acgtnACGTN\n]+"
        ocurrences = re.findall(FIND_FASTA_HEADER, data,re.MULTILINE)
        return ocurrences

    def read_byterange(self,ini,end):
        data = self.wsb.sb.read(ini,end)
        
    #Writes metadata from chunks to s3:
    #def write_metadata(self,obj,prev_header,index,chunk_index)
    def write_metadata(self,res,obj):
        self.storage.put_object(self.bucket_name, f'metadata/data{obj.part}','')
        for line in res:
            FIND_FASTA_HEADER = b"^>.*"
            pattern = re.compile(FIND_FASTA_HEADER,re.MULTILINE)
            ocurrences = re.findall(pattern, line)
            prev_data = self.storage.get_object(self.bucket_name, f'metadata/data{obj.part}',stream=True)
            if ocurrences != []:
                    self.storage.put_object(self.bucket_name, f'metadata/data{obj.part}',str(prev_data.read().decode("utf-8")) + str(ocurrences[0].decode("utf-8"))+" : "+str(len(line)) + '\n')
            else:
                    self.storage.put_object(self.bucket_name, f'metadata/data{obj.part}', str(len(line)) + '\n')

    #Waits for the previous chunk to be written, so it can read the last header from it.
    #def wait_for_prev_chunk(self,obj,prev_header,index,chunk_index)
    def get_chunk_waiting(self,mapper_num):

        storage = Storage()
        s3_client = storage.storage_handler.s3_client

        try:
            waiter = s3_client.get_waiter('object_exists')
            waiter.wait(Bucket=self.bucket_name, Key = f'fasta/data{mapper_num}',
            WaiterConfig={
                     'Delay': 5, 'MaxAttempts': 40})
            data  = self.storage.get_object(self.bucket_name, f'fasta/data{mapper_num}',stream=True)
        except Exception as e:
            print(f'fasta/data{mapper_num}' + " object doesn't exist")
            raise Exception( str(mapper_num)+":"+"Unexpected error in use_waiters_check_object_exists: " + e.__str__())
        
        return data

    #Writes the fasta chunk into S3
    #def write_fasta_chunk(self,res,obj,prev_header,index,chunk_index)
    def write_fasta_chunk(self,res,obj,prev_header,index,chunk_index,overlap):

        
        index = 0
        FASTA_HEADER = r">.+"
        prev_data = ""

        for data in res:
            data = data.decode("utf-8")
            found = re.search(FASTA_HEADER,data)
            if found != None:
                if obj.part == 0:
                    print("first chunk")
                    new_headers = re.sub(FASTA_HEADER,str(prev_header)+str(index)+"_"+str(chunk_index),data)
                    print(str(prev_header)+str(index)+"_"+str(chunk_index))
                    self.storage.put_object(self.bucket_name, f'fasta/data{obj.part}',str(prev_data)+str(new_headers))
                    prev_data = prev_data + new_headers
                    index += 1
                else:
                    print("not first chunk")
                    new_header = prev_header.split("_")
                    new_header[1] = str(int(new_header[1])+1)
                    new_header = "_".join(new_header)
                    new_headers = re.sub(FASTA_HEADER,str(new_header),data)
                    self.storage.put_object(self.bucket_name, f'fasta/data{obj.part}',str(prev_data)+str(new_headers))
                    prev_data = prev_data + new_headers 
                    prev_header = new_header
            else:
                print("part of another chunk")
                new_header = prev_header.split("_")
                new_header[2] = str(int(new_header[2])+1)
                new_header = "_".join(new_header)
                prev_data = str(new_header)+"\n"+str(data)
                self.storage.put_object(self.bucket_name, f'fasta/data{obj.part}',str(prev_data))
        

def print_chunk_info(obj):
    # obj is a CloudObject
    print(str("obj part :") + str(obj.part))



def run_worker(obj, storage,bucket_name,overlap):
    print_chunk_info(obj)
    #Mapper phase 1: Read data using lithops.utils.WrappedStreamingBodyPartition and chunk it into FASTA chunks based on 2 regex patterns.

    wsb = WrappedStreamingBodyPartitionFasta(obj.data_stream, storage, bucket_name)
    res = wsb.read()
    
    

    #Mapper phase 2: Read previous chunk in order to find the header last header
    #Prevent race condition between mappers, use of s3 waiters.
    FASTA_HEADER = r"^>.*"
    prev_mapper = obj.part - 1
    
    
    #First fasta chunk is always written, not caring about consistency.
    if obj.part == 0:
        wsb.write_fasta_chunk(res,obj,'>ch_',1,1,overlap)
    else:
    #Else, get the previous chunk header
        prev_chunk = wsb.get_chunk_waiting(prev_mapper).read().decode("utf-8")
        prev_header = re.findall(FASTA_HEADER, prev_chunk,re.MULTILINE)

        
    #Mapper phase 3: write the chunk into s3
        print("prev_header"+str(prev_header))
        wsb.write_fasta_chunk(res,obj,prev_header.pop(),1,1,overlap)

    return res


if __name__ == "__main__":
    fexec = lithops.FunctionExecutor(
        log_level='DEBUG',runtime='testcontainer', max_workers=1000, runtime_memory=4096)
    
    args = {'obj':'s3://sra-pub-src-2/SRR8774337/cleaned.fasta', 'bucket_name':'ayman-lithops-meta-cloudbutton-hutton', 'overlap':0}
    
    fexec.map(run_worker, args,
              obj_chunk_size=int(1*pow(2,20)))

    res = fexec.get_result()
    fexec.plot()
    fexec.job_summary()
 
    fexec.clean()

