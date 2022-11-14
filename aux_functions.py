import shutil
import os
from typing import List
from lithops import Storage
import re
import pickle
from lithops.utils import FuturesList
from varcall_arguments import Arguments

def create_fasta_chunk_for_runtime(storage: Storage, bucket: str, fasta: dict, byte_range: dict, folder: str, file_name: str):
    data = list(re.finditer(r">.+\n", storage.get_object(bucket=bucket, key=folder+file_name,
                extra_get_args={'Range': f"bytes={fasta['chunk'][0]['offset_head']}-{fasta['chunk'][0]['offset_base']}"}).decode('utf-8')))[0].group()
    base = storage.get_object(bucket=bucket, key=folder+file_name, extra_get_args=byte_range).decode('utf-8')
    if base[0:1] == '\n': # Data has already a '\n', (data == >...\n), avoid doble '\n'
        data += base[1::]
    else:
        data += base
    return data.encode('utf-8')

def copy_to_runtime(storage: Storage, bucket: str, folder: str, file_name: str, byte_range={}, fasta=None):
    """
    Copy file from S3 to local storage.
    """
    temp_file = "/tmp/" + file_name
    with open(temp_file, 'wb') as file:
        if fasta is not None: 
            file.write(create_fasta_chunk_for_runtime(storage, bucket, fasta, byte_range, folder, file_name))
        else:
            shutil.copyfileobj(storage.get_object(bucket=bucket, key=folder+file_name, stream=True, extra_get_args=byte_range), file)
    return temp_file


def copy_to_s3(storage: Storage, bucket: str, file_name: str, temp_to_s3: bool, folder_name: str = "") -> str:
    """
    Copy file from local storage to S3.
    """
    if temp_to_s3==True:
        destination_key= folder_name + os.path.basename(file_name)

        with open(file_name, 'rb') as file:
            data = file.read()
            storage.put_object(bucket, destination_key, body=data)
        return destination_key
    return ""


def load_cache(filename: str, args:Arguments) -> FuturesList:
    """
    Load a futures local file from previous execution
    """
    if os.path.isfile(f'/tmp/{args.execution_name}/{filename}') and args.checkpoints:
        file = open(f'/tmp/{args.execution_name}/{filename}', 'rb')
        futures = pickle.load(file)
        file.close()
        return futures
    else:
        return 0
    
def dump_cache(filename: str, futures: FuturesList, args:Arguments):
    """
    Store the lithops futures variable in local storage
    """
    if not os.path.exists(f'/tmp/{args.execution_name}'):
        os.makedirs(f'/tmp/{args.execution_name}')
    file = open(f'/tmp/{args.execution_name}/{filename}', 'wb')
    pickle.dump(futures, file)
    file.close()
    
def delete_files(storage: Storage, args:Arguments, cloud_prefixes: List[str] = [], local_files: List[str] = []):
    """
    Delete a list of cloud and local files
    """
    #Delete cloud files
    for prefix in cloud_prefixes:
        keys = storage.list_keys(args.bucket, prefix)
        for key in keys:
            storage.delete_object(args.bucket, key)
    
    #Delete local files        
    for file in local_files:
        if (os.path.isfile(file)):
            os.remove(file)