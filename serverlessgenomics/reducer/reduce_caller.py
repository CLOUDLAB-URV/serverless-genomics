from collections import defaultdict
from lithops import Storage
from typing import Tuple
import lithops as deflithops
import boto3

from ..parameters import PipelineRun, Lithops
from .reduce_functions import reduce_function
from ..stats import Stats
from ..utils import split_data_result

def keys_by_fasta_split(keys: Tuple[str]) -> dict:
    """
    Distribute the received keys by fasta split.

    Args:
        keys (Tuple[str]): List of keys

    Returns:
        dict: Dictionary where each key (0,1,2...) is the fasta split and the value is the list of keys
    """
    key_dict = defaultdict(list)
    for k in keys:
        fasta_split = k.split('/')[-1]
        fasta_split = fasta_split.split("fa")[-1]
        fasta_split = fasta_split[0]
        key_dict[str(fasta_split)].append(k)
    
    return key_dict


def create_multipart_keys(pipeline_params: PipelineRun) -> Tuple[str]:
    """
    Create the keys that will be used for the multipart uploads

    Args:
        pipeline_params (PipelineRun): Pipeline parameters

    Returns:
        Tuple[str]: List of keys
    """
    keys = []
    for i in range(pipeline_params.fasta_chunks):
        keys.append(f'tmp/{pipeline_params.run_id}/multipart_uploads/fa{i}.sinple')
    return keys


def create_multipart(pipeline_params: PipelineRun, key: str, storage: Storage) -> str:
    """
    Create a S3 multipart upload instance

    Args:
        pipeline_params (PipelineRun): Pipeline Parameters
        key (str): Key for the multipart upload
        storage (Storage): Lithops storage instance

    Returns:
        str: Multipart Upload ID
    """
    s3 = storage.get_client()
    mpu = s3.create_multipart_upload(
        Bucket=pipeline_params.storage_bucket,
        Key=key
    )
    return mpu['UploadId']


def distribute_indexes(pipeline_params: PipelineRun, keys: Tuple[str], storage: Storage) -> Tuple[Tuple[str]]:
    """
    Distribute the indexes between different reducers

    Args:
        pipeline_params (PipelineRun): Pipeline Parameters
        keys (Tuple[str]): Keys to the files generated in the map phase
        storage (Storage): Lithops storage instance

    Returns:
        Tuple[Tuple[str]]: Array where each position consists of the list of keys a reducer should take
    """
    mainStat, stat, subStat = Stats(), Stats(), Stats()
    mainStat.timer_start("call_distributeIndexes")
    s3 = storage.get_client()
    
    expression = "SELECT cast(s._2 as int) FROM s3object s"
    input_serialization = {'CSV': {'RecordDelimiter': '\n', 'FieldDelimiter': '\t'}, 'CompressionType': 'NONE'}
    
    count_indexes = {}
    
    # First we get the number of times each index appears
    stat.timer_start(expression)
    for key in keys:
        subStat.timer_start(key)
        resp = s3.select_object_content(
                    Bucket=pipeline_params.storage_bucket,
                    Key=key,
                    ExpressionType='SQL',
                    Expression=expression,
                    InputSerialization = input_serialization,
                    OutputSerialization = {'CSV': {}}
                )

        data = ""
        for event in resp['Payload']:
            if 'Records' in event:
                records = event['Records']['Payload'].decode("UTF-8")
                data = data + records
                
        data = data.split("\n")
        data.pop()  # Last value is empty
        
        int_indexes = list(map(int, data))

        for index in int_indexes:
            count_indexes[index] = count_indexes.get(index, 0) + 1
        subStat.timer_stop(key)
    stat.timer_stop(expression)
    stat.store_dictio(subStat.get_stats(), "subprocesses", expression)
    
    # Now we distribute the indexes depending on the max number of indexes we want each reducer to process
    MAX_INDEXES = 20_000_000
    workers_data = []
    indexes = 0
    
    for key in count_indexes:
        if indexes + count_indexes[key] < MAX_INDEXES:
            indexes += count_indexes[key]
            index = key
        else: # append the last index below max_index as end value in range, and start a new range.
            indexes = 0
            workers_data.append(index)
    workers_data.append(key)
    
    mainStat.timer_stop("call_distributeIndexes")
    mainStat.store_dictio(stat.get_stats(), "subprocesses", "call_distributeIndexes")
    
    return workers_data, mainStat.get_stats()

def create_iterdata_reducer(intermediate_keys: dict, distributed_indexes: Tuple[Tuple[str]], multipart_ids: Tuple[str], 
                            multipart_keys: Tuple[str], pipeline_params: PipelineRun) -> Tuple[dict]:
    """
    Create the iterdata for the reduce stage.

    Args:
        intermediate_keys (Tuple[str]): Keys distributed by fasta split
        distributed_indexes (Tuple[Tuple[str]]): Indexes that each reducer should process
        multipart_ids (Tuple[str]): Multipart Upload IDs
        multipart_keys (Tuple[str]): Multipart Upload Keys
        pipeline_params (PipelineRun): Pipeline Parameters

    Returns:
        Tuple[dict]: Reduce stage iterdata
    """
    iterdata = []

    for keys, indexes, mpu_id, mpu_key in zip(intermediate_keys, distributed_indexes, multipart_ids, multipart_keys):
        start = 1
        n_part = 1

        for index in indexes:
            data = {
                "keys" : intermediate_keys[keys], 
                "range" : { "start" : start,
                            "end" : int(index)
                        },
                "mpu_id" : mpu_id,
                "n_part" : n_part,
                "mpu_key" : mpu_key,
                "pipeline_params": pipeline_params
            }
            iterdata.append(data)
            n_part += 1
            start = int(index) + 1

    return iterdata


def final_merge(mpu_id: str, mpu_key: str, key: str, n_part: int, pipeline_params: PipelineRun, storage: Storage) -> dict:
    """
    Upload all the generated files by the reduce stage into one final file. This function will be mapped.

    Args:
        mpu_id (str): Multipart Upload ID
        mpu_key (str): Multipart Upload Key
        key (str): Key to the part file
        n_part (int): The number of the part file
        pipeline_params (PipelineRun): Pipeline Parameters
        storage (Storage): Lithops storage instance

    Returns:
        dict: Dictionary with the multipart upload settings
    """
    stat = Stats()
    stat.timer_start("call_finalMerge")

    sinple_out = storage.get_object(bucket=pipeline_params.storage_bucket, key=key)

    #Upload part
    s3 = storage.get_client()
    part = s3.upload_part(
        Body = sinple_out,
        Bucket = pipeline_params.storage_bucket,
        Key = mpu_key,
        UploadId = mpu_id,
        PartNumber = n_part
    )

    stat.timer_stop("call_finalMerge")

    return {"PartNumber" : n_part, "ETag" : part["ETag"], "mpu_id": mpu_id}, stat.get_stats()
    

def complete_multipart(keys: Tuple[str], mpu_ids: Tuple[str], parts: Tuple[dict], pipeline_params: PipelineRun, s3: Storage):
    """
    Complete a list of multipart uploads.

    Args:
        keys (Tuple[str]): Keys to the multipart uploads
        mpu_ids (Tuple[str]): IDs to the multipart uploads
        parts (Tuple[dict]): Parts of each multipart upload
        pipeline_params (PipelineRun): Pipeline Parameters
        s3 (Storage): Lithops storage instance
    """
    for key, mpu_id in zip(keys, mpu_ids):
        mpu_part = []
        remove = 0

        for part in parts:
            if mpu_id == part['mpu_id']:
                mpu_part.append({"PartNumber" : part["PartNumber"], "ETag" : part["ETag"]})
                remove = remove + 1
            else:
                break

        s3.complete_multipart_upload(
            Bucket = pipeline_params.storage_bucket,
            Key = key,
            UploadId = mpu_id,
            MultipartUpload = {"Parts": mpu_part}
        )

        #for _ in range(remove):
        #    parts.pop(0)
        parts = parts[remove:]

def finish(key: str, mpu_id: str, parts: Tuple[dict], pipeline_params: PipelineRun, s3: Storage):
    """
    Complete the final multipart upload

    Args:
        key (str): Multipart upload Key
        mpu_id (str): Multipart upload ID
        parts (Tuple[dict]): Multipart upload parts
        pipeline_params (PipelineRun): Pipeline parameters
        s3 (Storage): Lithops storage instance
    """
    mpu_part = []
    remove = 0

    for part in parts:
        if mpu_id == part['mpu_id']:
            mpu_part.append({"PartNumber" : part["PartNumber"], "ETag" : part["ETag"]})
            remove = remove + 1
        else:
            break

    s3.complete_multipart_upload(
        Bucket = pipeline_params.storage_bucket,
        Key = key,
        UploadId = mpu_id,
        MultipartUpload = {"Parts": mpu_part}
    )

def run_reducer(pipeline_params: PipelineRun, lithops: Lithops, mapper_output):
    subStat = Stats()
    # 1 Organize the keys generated by the map phase by fasta split
    intermediate_keys = keys_by_fasta_split(mapper_output)

    # 2 Create the keys for the multipart uploads
    multipart_keys = create_multipart_keys(pipeline_params)
 
    # 3 Create the multipart uploads and get their IDs
    multipart_ids = []
    for key in multipart_keys:
        multipart_ids.append(create_multipart(pipeline_params, key, lithops.storage))

    # 4 Select the indexes that each reducer will process
    indexes_iterdata = []
    for fasta in intermediate_keys:
        data = {
            'pipeline_params': pipeline_params,
            'keys': intermediate_keys[fasta]
        }
        indexes_iterdata.append(data)
    
    subStat.timer_start('distribute_indexes')
    distributed_indexes = lithops.invoker.map(distribute_indexes, indexes_iterdata)
    subStat.timer_stop('distribute_indexes')
    distributed_indexes, timers = split_data_result(distributed_indexes)
    subStat.store_dictio(timers, "subprocesses", "distribute_indexes") 
    
    # 5 Launch the reducers
    reducer_iterdata = create_iterdata_reducer(intermediate_keys, distributed_indexes, multipart_ids, multipart_keys, pipeline_params)
    subStat.timer_start('reduce_function')
    reducer_output = lithops.invoker.map(reduce_function, reducer_iterdata)
    subStat.timer_stop('reduce_function')
    reducer_output, timers = split_data_result(reducer_output)
    subStat.store_dictio(timers, "subprocesses", "reduce_function")

   
    
    # 6 Complete the multipart uploads that the reducers created
    complete_multipart(multipart_keys, multipart_ids, reducer_output, pipeline_params, lithops.storage.storage_handler.s3_client)
    
    # 7 Create a multipart upload key and ID for the final file
    final_sinple_key = f'tmp/{pipeline_params.run_id}/final.alignment'
    final_id = create_multipart(pipeline_params, final_sinple_key, lithops.storage)
    
    # 8 Merge files created in stage 6 into one single file
    n_parts = len(multipart_keys)
    part = 1
    merge_iterdata = []
    while part <= n_parts:
        data = {
            "mpu_id": final_id,
            "mpu_key": final_sinple_key,
            "key": multipart_keys[part-1],
            "n_part": part,
            "pipeline_params": pipeline_params
        }
        merge_iterdata.append(data)
        part += 1
    
    subStat.timer_start('final_merge')
    final_merge_results = lithops.invoker.map(final_merge, merge_iterdata)
    subStat.timer_stop('final_merge')
    final_merge_results, timers = split_data_result(final_merge_results)
    subStat.store_dictio(timers, "subprocesses", "final_merge")
       

    # 8 Complete the previous multipart upload
    finish(final_sinple_key, final_id, final_merge_results, pipeline_params, lithops.storage.storage_handler.s3_client)

    return subStat