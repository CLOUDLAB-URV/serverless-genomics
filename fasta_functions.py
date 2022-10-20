import os
from lithops import Storage
import subprocess as sp

def prepare_fasta(cloud_adr, runtime, BUCKET_NAME, fasta_folder, fasta_file, split_fasta_folder, fasta_chunk_size, fasta_chunks_prefix, fasta_line_overlap):
    """
    Try to find fasta chunks in storage. If they are not found proceed to partition the orignal fasta file.
    TODO Replace partitioner with the new implementation.
    """
    
    storage = Storage()
    fasta_list = []
    
    # Check if chunks exist in storage
    try:
        fasta_list = storage.list_keys(BUCKET_NAME, prefix=split_fasta_folder)
    except:
        print("split fasta folder empty / not found")


    # If chunks don't exist in storage, split fasta file.
    if not fasta_list or fasta_list == []:
        print("splitting fasta file and storing chunks")
        split_fastafile(runtime, BUCKET_NAME, fasta_folder, fasta_file, split_fasta_folder, fasta_chunk_size, fasta_line_overlap)
        try:
            print("Finding newly created fasta chunks")
            print("Searching " + BUCKET_NAME + "/" + split_fasta_folder + fasta_chunks_prefix)
            print("Finding in storage: \n" + str(storage.list_keys(BUCKET_NAME, prefix=split_fasta_folder + fasta_chunks_prefix)))
            fasta_list = storage.list_keys(BUCKET_NAME, prefix=split_fasta_folder + fasta_chunks_prefix)
        except:
            print("error generating fasta chunks")
    return fasta_list


def split_fastafile(runtime, bucket_name, fasta_folder, fasta_file, split_fasta_folder, fasta_chunk_size, chunk_overlap):
    """
    Splits fasta file into chunks depending on the desired chunk size
    TODO Replace partitioner with the new implementation.
    """

    args = "python FastaPartitioner.py --data_location "+"s3"+"://"+bucket_name+"/"+fasta_folder+fasta_file+" --chunk_size "+str(fasta_chunk_size)+" --overlap "+str(chunk_overlap)+" --mybucket "+bucket_name+" --fasta_folder "+split_fasta_folder+" --runtime "+runtime
    args = args.split(" ")

    total_objects = sp.check_output(args).decode('utf-8')

    i = 0

    storage = Storage()
    data = ""
    chr_title = ""
    ini = 1
    while i <= int(total_objects):

        ret = storage.get_object(bucket_name, f'cache/obj{i}.data', stream = True).read().decode('utf-8')
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
    storage.put_object(bucket_name, f'cache/hashtable.data',data)
    print("end of *split_fasta* function")