import os
from typing import List
import pickle
from lithops.utils import FuturesList
from lithops import Storage
from .parameters import PipelineRun


def copy_to_s3(storage: Storage, bucket: str, file_name: str, temp_to_s3: bool, folder_name: str = "") -> str:
    """
    Copy file from local storage to S3.
    """
    if temp_to_s3 == True:
        destination_key = folder_name + os.path.basename(file_name)

        with open(file_name, 'rb') as file:
            data = file.read()
            storage.put_object(bucket, destination_key, body=data)
        return destination_key
    return ""


def load_cache(filename: str, args: PipelineRun) -> FuturesList:
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


def dump_cache(filename: str, futures: FuturesList, args: PipelineRun):
    """
    Store the lithops futures variable in local storage
    """
    if not os.path.exists(f'/tmp/{args.execution_name}'):
        os.makedirs(f'/tmp/{args.execution_name}')
    file = open(f'/tmp/{args.execution_name}/{filename}', 'wb')
    pickle.dump(futures, file)
    file.close()


def delete_files(storage: Storage, args: PipelineRun, cloud_prefixes: List[str] = [], local_files: List[str] = []):
    """
    Delete a list of cloud and local files
    """
    # Delete cloud files
    for prefix in cloud_prefixes:
        keys = storage.list_keys(args.storage_bucket, prefix)
        for key in keys:
            storage.delete_object(args.storage_bucket, key)

    # Delete local files
    for file in local_files:
        if os.path.isfile(file):
            os.remove(file)
