from lithops import Storage
import subprocess as sp
from varcall_arguments import Arguments

def prepare_fasta(args: Arguments, fasta_chunks_prefix: str):
    """
    Try to find fasta chunks in storage. If they are not found proceed to partition the orignal fasta file.
    TODO Replace partitioner with the new implementation.
    """
    
    storage = Storage()
    fasta_list = []
    
    # Check if chunks exist in storage
    try:
        fasta_list = storage.list_keys(args.bucket, prefix=args.split_fasta_folder)
    except:
        print("split fasta folder empty / not found")


    # If chunks don't exist in storage, split fasta file.
    if not fasta_list or fasta_list == []:
        print("splitting fasta file and storing chunks")
        split_fastafile(args.runtime_id, args.bucket, args.fasta_folder, args.fasta_file, args.split_fasta_folder, args.fasta_chunk_size, args.fasta_char_overlap)
        try:
            print("Finding newly created fasta chunks")
            print("Searching " + args.bucket + "/" + args.split_fasta_folder + fasta_chunks_prefix)
            print("Finding in storage: \n" + str(storage.list_keys(args.bucket, prefix=args.split_fasta_folder + fasta_chunks_prefix)))
            fasta_list = storage.list_keys(args.bucket, prefix=args.split_fasta_folder + fasta_chunks_prefix)
        except:
            print("error generating fasta chunks")
    return fasta_list

# runtime, bucket_name, fasta_folder, fasta_file, split_fasta_folder, fasta_chunk_size, chunk_overlap
def split_fastafile(args: Arguments):
    """
    Splits fasta file into chunks depending on the desired chunk size
    TODO Replace partitioner with the new implementation.
    """

    split_args = "python FastaPartitioner.py --data_location "+"s3"+"://"+args.bucket+"/"+args.fasta_folder+args.fasta_file+" --chunk_size "+str(args.fasta_chunk_size)+" --overlap "+str(args.fasta_char_overlap)+" --mybucket "+args.bucket+" --fasta_folder "+args.split_fasta_folder+" --runtime "+args.runtime_id
    split_args = split_args.split(" ")

    total_objects = sp.check_output(split_args).decode('utf-8')

    i = 0

    storage = Storage()
    data = ""
    chr_title = ""
    ini = 1
    while i <= int(total_objects):

        ret = storage.get_object(args.bucket, f'cache/obj{i}.data', stream = True).read().decode('utf-8')
        for line in ret.splitlines():
            if line == '\n' or line == '':
                pass
            elif '-' in line and ',' in line:
                d = line.split(',')
                range = d[1].split('-')
                abs_pos = int(range[1]) -  int(range[0])+1
                range = '\t'.join(range)
                data =   data + d[0].replace('>','')+ '\t' +chr_title[1].replace('>','') + '\t' + str(ini) + '\t' +  str(abs_pos)  + '\n'
                ini = abs_pos
            elif '-' in line:
                range = line.split('-')
                abs_pos = int(range[1]) -  int(range[0])+1
                data =  data + chr_title[0].replace('>','') + '\t' + chr_title[1].replace('>','') + '\t' + str(ini) + '\t' +  str(abs_pos) +'\n'
                ini = abs_pos
            else:
                chr_title = line.split(',')
                ini = 1
        i= i+1

    print(data)
    storage.put_object(args.bucket, f'cache/hashtable.data',data)
    print("end of *split_fasta* function")