from collections import defaultdict
import os
from subprocess import Popen, PIPE, STDOUT
from typing import Tuple
from time import time
from sys import getsizeof
from pprint import pprint

from lithops import Storage
from ..pipeline import PipelineParameters, PipelineRun
from ..stats import Stats


def reduce_function(keys, range, mpu_id, n_part, mpu_key, pipeline_params: PipelineParameters, storage: Storage):
    stats = Stats()
    stats.start_timer("function")
    stats.set_value("keys", keys)
    stats.set_value("range", range)
    stats.set_value("mpu_id", mpu_id)
    stats.set_value("n_part", n_part)
    stats.set_value("mpu_key", mpu_key)

    s3 = storage.get_client()

    # Change working directory to /tmp
    wd = os.getcwd()
    os.chdir("/tmp")

    # S3 SELECT query to get the rows where the second column is in the selected range
    expression = "SELECT * FROM s3object s WHERE cast(s._2 as int) BETWEEN %s AND %s" % (range["start"], range["end"])
    input_serialization = {
        "CSV": {"RecordDelimiter": "\n", "FieldDelimiter": "\t"},
        "CompressionType": "NONE",
    }

    # Execute S3 SELECT
    mpileup_data = ""
    for k in keys:
        key_stat = Stats()
        key_stat.set_value("key", k)
        try:
            with key_stat.timeit("s3_select"):
                resp = s3.select_object_content(
                    Bucket=pipeline_params.storage_bucket,
                    Key=k,
                    ExpressionType="SQL",
                    Expression=expression,
                    InputSerialization=input_serialization,
                    OutputSerialization={"CSV": {"FieldDelimiter": "\t"}},
                )
        except:
            raise ValueError("ERROR IN KEY: " + k)

        for event in resp["Payload"]:
            if "Records" in event:
                records = event["Records"]["Payload"].decode("UTF-8")
                key_stat.set_value("data_size", len(records))
                mpileup_data = mpileup_data + records
                del records

        stats.set_value(k, key_stat.dump_dict())

    stats.set_value("mpileup_data_size", len(mpileup_data))

    # Str -> bytes
    mpileup_data = mpileup_data.encode("UTF-8")

    # Execute the script to merge and reduce
    with stats.timeit("mpileup_merge_reduce"):
        p = Popen(
            ["bash", "/function/bin/mpileup_merge_reducev3.sh", "/function/bin/", "75%"],
            stdout=PIPE,
            stdin=PIPE,
            stderr=PIPE,
        )
        sinple_out = p.communicate(input=mpileup_data)[0]
    sinple_out = sinple_out.decode("UTF-8")

    # Upload part
    with stats.timeit("upload_part"):
        part = s3.upload_part(
            Body=sinple_out, Bucket=pipeline_params.storage_bucket, Key=mpu_key, UploadId=mpu_id, PartNumber=n_part
        )

    os.chdir(wd)
    stats.stop_timer("function")
    return {"PartNumber": n_part, "ETag": part["ETag"], "mpu_id": mpu_id}, stats


def distribute_indexes(
    pipeline_params: PipelineParameters, fasta_chunk: int, keys: Tuple[str], storage: Storage
) -> Tuple[Tuple[str]]:
    """
    Distribute the indexes between different reducers

    Args:
        pipeline_params (PipelineParameters): Pipeline Parameters
        keys (Tuple[str]): Keys to the files generated in the map phase
        storage (Storage): Lithops storage instance

    Returns:
        Tuple[Tuple[str]]: Array where each position consists of the list of keys a reducer should take
    """
    stats = Stats()
    stats.start_timer("function")
    stats.set_value("fasta_chunk", fasta_chunk)
    stats.set_value("keys", keys)

    s3 = storage.get_client()

    expression = "SELECT cast(s._2 as int) FROM s3object s"
    input_serialization = {
        "CSV": {"RecordDelimiter": "\n", "FieldDelimiter": "\t"},
        "CompressionType": "NONE",
    }

    count_indexes = {}

    # First we get the number of times each index appears
    with stats.timeit("s3_select"):
        for key in keys:
            key_stats = Stats()
            key_stats.set_value("key", key)
            with key_stats.timeit("s3_select"):
                resp = s3.select_object_content(
                    Bucket=pipeline_params.storage_bucket,
                    Key=key,
                    ExpressionType="SQL",
                    Expression=expression,
                    InputSerialization=input_serialization,
                    OutputSerialization={"CSV": {}},
                )

            data = ""
            for event in resp["Payload"]:
                if "Records" in event:
                    records = event["Records"]["Payload"].decode("UTF-8")
                    data = data + records
            data = data.split("\n")
            data.pop()  # Last value is empty

            key_stats.set_value("data_size", len(data))

            int_indexes = list(map(int, data))

            for index in int_indexes:
                count_indexes[index] = count_indexes.get(index, 0) + 1

            stats.set_value(key, key_stats.dump_dict())

    # Now we distribute the indexes depending on the max number of indexes we want each reducer to process
    with stats.timeit("distribute_indexes"):
        MAX_INDEXES = 20_000_000
        workers_data = []
        indexes = 0

        for key in count_indexes:
            if indexes + count_indexes[key] < MAX_INDEXES:
                indexes += count_indexes[key]
                index = key
            else:  # append the last index below max_index as end value in range, and start a new range.
                indexes = 0
                workers_data.append(index)
        workers_data.append(key)

    stats.stop_timer("function")
    return workers_data, stats


def final_merge(
    mpu_id: str,
    mpu_key: str,
    key: str,
    n_part: int,
    pipeline_params: PipelineParameters,
    storage: Storage,
) -> dict:
    """
    Upload all the generated files by the reduce stage into one final file. This function will be mapped.

    Args:
        mpu_id (str): Multipart Upload ID
        mpu_key (str): Multipart Upload Key
        key (str): Key to the part file
        n_part (int): The number of the part file
        pipeline_params (PipelineParameters): Pipeline Parameters
        storage (Storage): Lithops storage instance

    Returns:
        dict: Dictionary with the multipart upload settings
    """
    stats = Stats()
    stats.start_timer("function")
    stats.set_value("mpu_id", mpu_id)
    stats.set_value("mpu_key", mpu_key)
    stats.set_value("key", key)
    stats.set_value("n_part", n_part)

    with stats.timeit("download_sinple_out"):
        sinple_out = storage.download_file(bucket=pipeline_params.storage_bucket, key=key)

    stats.set_value("sinple_out_size", len(sinple_out))

    # Upload part
    s3 = storage.get_client()
    with stats.timeit("upload_part"):
        part = s3.upload_part(
            Body=sinple_out, Bucket=pipeline_params.storage_bucket, Key=mpu_key, UploadId=mpu_id, PartNumber=n_part
        )

    stats.stop_timer("function")
    return {"PartNumber": n_part, "ETag": part["ETag"], "mpu_id": mpu_id}, stats


def finish(
    key: str,
    mpu_id: str,
    parts: Tuple[dict],
    pipeline_params: PipelineParameters,
    s3: Storage,
):
    """
    Complete the final multipart upload

    Args:
        key (str): Multipart upload Key
        mpu_id (str): Multipart upload ID
        parts (Tuple[dict]): Multipart upload parts
        pipeline_params (PipelineParameters): Pipeline parameters
        s3 (Storage): Lithops storage instance
    """
    mpu_part = []
    remove = 0

    for part in parts:
        if mpu_id == part["mpu_id"]:
            mpu_part.append({"PartNumber": part["PartNumber"], "ETag": part["ETag"]})
            remove = remove + 1
        else:
            break

    s3.complete_multipart_upload(
        Bucket=pipeline_params.storage_bucket,
        Key=key,
        UploadId=mpu_id,
        MultipartUpload={"Parts": mpu_part},
    )


def complete_multipart(
    keys: Tuple[str],
    mpu_ids: Tuple[str],
    parts: Tuple[dict],
    pipeline_params: PipelineParameters,
    s3: Storage,
):
    """
    Complete a list of multipart uploads.

    Args:
        keys (Tuple[str]): Keys to the multipart uploads
        mpu_ids (Tuple[str]): IDs to the multipart uploads
        parts (Tuple[dict]): Parts of each multipart upload
        pipeline_params (PipelineParameters): Pipeline Parameters
        s3 (Storage): Lithops storage instance
    """
    for key, mpu_id in zip(keys, mpu_ids):
        mpu_part = []
        remove = 0

        for part in parts:
            if mpu_id == part["mpu_id"]:
                mpu_part.append({"PartNumber": part["PartNumber"], "ETag": part["ETag"]})
                remove = remove + 1
            else:
                break

        s3.complete_multipart_upload(
            Bucket=pipeline_params.storage_bucket,
            Key=key,
            UploadId=mpu_id,
            MultipartUpload={"Parts": mpu_part},
        )

        parts = parts[remove:]


def keys_by_fasta_split(keys: Tuple[str]) -> dict:
    """
    Distribute the received keys by fasta split.

    Args:
        keys (Tuple[str]): List of keys

    Returns:
        dict: Dictionary where each key (0,1,2...) is the fasta split and the value is the list of keys
    """
    key_dict = defaultdict(list)
    for key in keys:
        worker_id = key.split("/")[-2]
        fasta_chunk_id = worker_id.split("-")[0]
        fasta_chunk_num = int(fasta_chunk_id.replace("fa", ""))
        key_dict[fasta_chunk_num].append(key)

    return key_dict


def create_multipart_keys(pipeline_params: PipelineParameters, pipeline_run: PipelineRun) -> Tuple[str]:
    """
    Create the keys that will be used for the multipart uploads

    Args:
        pipeline_params (PipelineParameters): Pipeline parameters

    Returns:
        Tuple[str]: List of keys
    """
    keys = []
    for i in range(pipeline_params.fasta_chunks):
        if (pipeline_params.fasta_chunk_range is None) or (i in pipeline_params.fasta_chunk_range):
            keys.append(f"tmp/{pipeline_run.run_id}/multipart_uploads/fa{i}.sinple")
    return keys


def create_multipart(pipeline_params: PipelineParameters, key: str, storage: Storage) -> str:
    """
    Create a S3 multipart upload instance

    Args:
        pipeline_params (PipelineParameters): Pipeline Parameters
        key (str): Key for the multipart upload
        storage (Storage): Lithops storage instance

    Returns:
        str: Multipart Upload ID
    """
    s3 = storage.get_client()
    mpu = s3.create_multipart_upload(Bucket=pipeline_params.storage_bucket, Key=key)
    return mpu["UploadId"]
