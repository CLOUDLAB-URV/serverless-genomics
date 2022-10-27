from lithops import Storage
import subprocess as sp
from varcall_arguments import Arguments
import pathlib

def prepare_fasta(args: Arguments):
    """
    Try to find fasta chunks in storage. If they are not found proceed to partition the orignal fasta file.
    """
    
    storage = Storage()
    fasta_index = args.fasta_folder_index + pathlib.Path(args.fasta_file).stem + f'_{args.fasta_workers}.fai'   
    try:
        print("Searching " + args.bucket + "/" + fasta_index)
        list_files_index = storage.list_keys(args.bucket, prefix=args.fasta_workers)
        print("Finding in storage: \n" + str(list_files_index))
        storage.head_object(args.bucket, fasta_index)
        
    except:
        print("index fasta folder empty / not found, then generating index fasta file and storing it")
        create_fastafile_index(args.bucket, args.fasta_folder, args.fasta_file, args.fasta_folder_index, args.fasta_workers)
        try:
            print("Searching " + args.bucket + "/" + args.fasta_folder_index)
            print("Finding in storage: \n" + str(storage.list_keys(args.bucket, prefix=args.fasta_folder_index)))
            storage.head_object(args.bucket, fasta_index)
        except:
            print("Error generating fasta file index")
    return fasta_index

def create_fastafile_index(bucket_name, fasta_folder, fasta_file, fasta_folder_index, fasta_n_workers):
    """
    Splits fasta file into chunks depending on the desired chunk size
    """
    storage = Storage()

    print("start of *create_fastafile_index* function")
    print("--> Generating fasta file index")

    # prints all files
    args = "python fastaPartitionerIndex.py --key "+fasta_folder+fasta_file+" --workers "+str(fasta_n_workers)+" --mybucket "+bucket_name+" --fasta_folder "+fasta_folder_index

    args = args.split(" ")
    # Modified, total_objects = sp.check_output(args).decode('utf-8')
    sp.check_output(args)

    print("end of *create_fastafile_index* function")
