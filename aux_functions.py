import shutil
import os
from lithops import Storage

def copy_to_runtime(storage: Storage, bucket: str, folder: str, file_name: str, byte_range=None):
    """
    Copy file from S3 to local storage.
    """
    extra_get_args = {'Range': f'bytes={byte_range}'} if byte_range else {}
    obj_stream = storage.get_object(bucket=bucket, key=folder+file_name, stream=True, extra_get_args=extra_get_args)
    temp_file = "/tmp/" + file_name
    with open(temp_file, 'wb') as file:
        shutil.copyfileobj(obj_stream, file)
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