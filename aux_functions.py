import shutil
import os
import re
import time
import subprocess as sp
import pathlib
import tempfile

CURRENT_PATH = str(pathlib.Path(__file__).parent.resolve())
LOCAL_TMP = os.path.realpath(tempfile.gettempdir())


def copy_to_runtime(storage, bucket, folder, file_name, byte_range=None):
    """
    Copy file from S3 to local storage.
    """
    extra_get_args = {'Range': f'bytes={byte_range}'} if byte_range else {}
    obj_stream = storage.get_object(bucket=bucket, key=folder+file_name, stream=True, extra_get_args=extra_get_args)
    temp_file = "/tmp/" + file_name
    with open(temp_file, 'wb') as file:
        shutil.copyfileobj(obj_stream, file)
    return temp_file


def copy_to_s3(storage, BUCKET_NAME, file_name, temp_to_s3, folder_name=""):
    """
    Copy file from local storage to S3.
    """
    if temp_to_s3==True:
        destination_key= folder_name + os.path.basename(file_name)

        with open(file_name, 'rb') as file:
            data = file.read()
            storage.put_object(BUCKET_NAME, destination_key, body=data)
        return destination_key
    return 0