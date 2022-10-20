
import chunk
import re
import lithops
import argparse
import sys
from PriceEstimator import PriceEstimator

#### Aded by Sara
BUCKET_NAME = "genomics-sara"
NAME_FASTA = 'TriTrypDB-9.0_TbruceiTREU927_Genome_MbChr' # Added by Sara
ENDPOINT = 'http://192.168.2.3:9000/'
WORKERS = 1000000
OBJ_LOCATION = f'{ENDPOINT}{BUCKET_NAME}/fasta/{NAME_FASTA}.fasta'
OBJ_KEY = f'fasta/{NAME_FASTA}.fasta'
####

parser = argparse.ArgumentParser(
    description='Fasta partitioner, version 4.2. Takes a fasta file and chunks it into multiple parts.',
)

parser.add_argument('--mybucket',
                    help='Your bucket name, where you want data to be located, ex: ayman-lithops-meta-cloudbutton-hutton')
parser.add_argument('--chunk_size',
                    help='Size of the .fasta you want in bytes')
parser.add_argument('--overlap',
                    help='Overlap of the fasta chunks')
parser.add_argument('--data_location',
                    help=f'Location of the fasta file, ex: s3://ayman-lithops-meta-cloudbutton-hutton/fasta/hg19.fa')
parser.add_argument('--runtime',
                    help='Runtime to use for FunctionExecutor')
parser.add_argument('--fasta_folder',
                    help='Fasta folder within your bucket to store the chunks, ex: fasta/')
args = parser.parse_args()


obj = args.data_location
my_bucket_name = args.mybucket
worker_chunk_size = int(args.chunk_size)
ovlp = int(args.overlap)
fasta_folder = args.fasta_folder
runtime = args.runtime





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
                    content =  content + prev +',' + str(start-1)+ '\n'
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


        self.storage.put_object(self.my_bucket_name, f'cache/obj{obj.part-1}.data',content)

        return f'cache/obj{obj.part}.data'


    def generate_index_file_titles(self,total_obj_size):
        i=0
        sub_index = ""
        chr_index = 0
        chr_subindex = 0
        obj_key = OBJ_KEY

        total_chunks = int(total_obj_size)
    
        while i < total_chunks+1:
            data = self.storage.get_object(self.my_bucket_name, f'cache/obj{i}.data', stream=True)
            data_lines = data.read().decode('utf-8').splitlines()
            num_lines = len(data_lines)
            title = False
            for el in data_lines:
                if 'x' in el:
                    start = int(el.split(',')[0])
                    end = int(el.split(',')[1].replace('x',''))
                    chunk = self.storage.get_object(self.my_bucket_name, obj_key,stream=True,extra_get_args = {'Range': f'bytes={start}-{int(end)-1}'}).read().decode('utf-8')
                    chr_index = chr_index + 1
                    chr_subindex = 1
                    sub_index = sub_index + ('>'+str(chr_index) + '_' + str(chr_subindex) +','+ chunk ) + '\n'
                    title = True
                else:
                    if title == False:
                        chr_subindex = chr_subindex + 1
                        start = el.split(',')[0]
                        end = el.split(',')[1]
                        sub_index = sub_index + ('>'+str(chr_index) + '_' + str(chr_subindex) +','+ str(start) + '-' + str(end)) + '\n'
                    else:
                        start = el.split(',')[0]
                        end = el.split(',')[1]
                        sub_index = sub_index + (str(start) + '-' + str(end)) + '\n'
                        title = False

            self.storage.put_object(self.my_bucket_name, f'cache/obj{i}.data',str(sub_index))
            sub_index = ""
            i= i+1


def print_chunk_info(obj):
    # obj is a CloudObject
    print("===================")
    print(obj)
    print(obj.key)
    print('part:'+str(obj.part))
    print(obj.data_byte_range)
    print(obj.chunk_size)  # in bytes
    print("===================")

def run_worker_metadata(obj,storage,my_bucket_name,total_obj_size,obj_size,overlap):
    #print_chunk_info(obj)
    partitioner = PartitionFasta(storage,my_bucket_name,obj)
    chunks = partitioner.generate_chunks(obj,total_obj_size,obj_size,overlap)
    return chunks

def generate_index_file(obj,storage,my_bucket_name,total_obj_size):
    partitioner = PartitionFasta(storage,my_bucket_name,obj)
    partitioner.generate_index_file_titles(total_obj_size)

def generate_chunks_corrected_index(id,r,storage):

    location = obj.split('//')[1]
    obj_bucket = BUCKET_NAME
    obj_key = location.split('/')[1:]
    obj_key = OBJ_KEY # Change-me
    obj_name = NAME_FASTA
    
    chunk_name = r.split('/')[-1]
    chunk_name = chunk_name.split('.')[0]
    chunk_num = int(chunk_name.split('j')[1])-1
    chunk_name = 'cache/obj' + str(chunk_num) + '.data'
    r = chunk_name
    
    start=0
    end=0
    ret_data = storage.get_object(my_bucket_name, r, stream=True).read().decode('utf-8')
    data = ""

    for el in ret_data.splitlines():
        if ',' in el and '-' in el:
            data = data + str(el.split(',')[0]) + '\n'
            start = int(el.split(',')[1].split('-')[0])
            end =  int(el.split(',')[1].split('-')[1])
            data = data + str(storage.get_object(obj_bucket, obj_key, stream=True,extra_get_args = {'Range': f'bytes={start}-{end}'}).read().decode('utf-8')) + '\n'
        elif ',' in el:
            data = data + str(el.split(',')[0]) + '\n'
        elif '-' in el:
            start = int(el.split('-')[0])
            end = int(el.split('-')[1])
            data = data + str(storage.get_object(obj_bucket, obj_key, stream=True,extra_get_args = {'Range': f'bytes={start}-{end}'}).read().decode('utf-8')) + '\n'
    storage.put_object(my_bucket_name, f'{fasta_folder}{obj_name}_{worker_chunk_size}split_{id}.fasta',data)
    data = ""






#s3://sra-pub-src-2/SRR8774337/cleaned.fasta

if __name__ == "__main__":
    fexec = lithops.FunctionExecutor(log_level='DEBUG', runtime_memory=16384)

    
    obj = OBJ_LOCATION  # Change-me
    
    my_bucket_name = BUCKET_NAME # Change-me
    fasta_key = OBJ_KEY # Change-me
    worker_chunk_size = WORKERS
    ovlp = 300
    
    arr = []
    location = obj.split('//')[1].split('/')
    obj_location = location[1]
    for_head = location[1:]
    for_head = fasta_key
    data_bucket_name = my_bucket_name
    storage = lithops.Storage()

    fasta = storage.head_object(data_bucket_name, for_head)
    total_fasta_size = int(fasta['content-length'])/worker_chunk_size
    seq_name = location[-1].split('.')[0]

    fexec.map(run_worker_metadata, {'my_bucket_name':my_bucket_name,'obj':obj,'total_obj_size':int(total_fasta_size),'obj_size':fasta['content-length'], 'overlap':ovlp}, obj_chunk_size=worker_chunk_size)
    iterdata = fexec.get_result()
    arr.append([len(fexec.futures),4])

    #First mapper doesn't always need to be executed (Caching feature)
    fexec.call_async(generate_index_file,{'my_bucket_name':my_bucket_name,'obj':obj,'total_obj_size':int(total_fasta_size)},timeout=60)
    fexec.wait()
    arr.append([ len(fexec.futures),3])



    fexec.map(generate_chunks_corrected_index, iterdata)
    fexec.wait()
    arr.append([len(fexec.futures),5])

    # pEst = PriceEstimator(fexec)
    #print('price lambda:'+ str(pEst.lambda_calc(4)))

    # pEst.arr = arr
    #print(pEst.s3_calc())


    sys.stdout.write(str(int(total_fasta_size)))
