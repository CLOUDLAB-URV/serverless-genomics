from lithops import Storage
import subprocess as sp
from serverlessgenomics.parameters import PipelineParameters
import pathlib

from .fasta_partitioner_index import FastaPartitionerCaller
import re

def prepare_fasta(args: PipelineParameters):
    """
    Try to find fasta chunks in storage. If they are not found proceed to partition the orignal fasta file.
    """

    storage = Storage()
    fasta_index = args.fasta_folder_index + pathlib.Path(args.fasta_file).stem + '.fai'
    try:
        print("Searching " + args.bucket + "/" + fasta_index)
        storage.head_object(args.bucket, fasta_index)

    except:
        print("Index fasta folder empty / not found, then generating index fasta file and storing it")
        fasta_partitioner = FastaPartitionerCaller(args.bucket)
        fasta_partitioner(args.fasta_folder + args.fasta_file, int(args.fasta_workers), args.fasta_folder_index)
        try:
            print("Searching " + args.bucket + "/" + args.fasta_folder_index)
            print("Finding in storage: \n" + str(storage.list_keys(args.bucket, prefix=args.fasta_folder_index)))
            storage.head_object(args.bucket, fasta_index)
        except:
            print("Error generating fasta file index")
    return fasta_index


def create_fasta_chunk_for_runtime(storage: Storage, bucket: str, fasta: dict, byte_range: dict, folder: str, file_name: str):
    data = list(re.finditer(r">.+\n", storage.get_object(bucket=bucket, key=folder+file_name,
                extra_get_args={'Range': f"bytes={fasta['chunk']['offset_head']}-{fasta['chunk']['offset_base']}"}).decode('utf-8')))[0].group()
    base = storage.get_object(bucket=bucket, key=folder+file_name, extra_get_args=byte_range).decode('utf-8')
    if base[0:1] == '\n': # Data has already a '\n', (data == >...\n), avoid doble '\n'
        data += base[1::]
    else:
        data += base
    return data.encode('utf-8')

    