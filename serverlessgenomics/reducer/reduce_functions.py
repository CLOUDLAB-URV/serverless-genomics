import os
import subprocess as sp

from lithops import Storage

from ..parameters import PipelineRun, Lithops
from ..stats import Stats

def reduce_function(keys, range, mpu_id, n_part, mpu_key, pipeline_params: PipelineRun, storage: Storage):
    mainStat, stat, subStat = Stats(), Stats(), Stats()
    mainStat.timer_start("call_reduceFunction")
    s3 = storage.get_client()
     
    # Change working directory to /tmp
    wd = os.getcwd()
    os.chdir("/tmp")
    
    # File where we will store the data we get from the SELECT queries
    temp_mpileup = '/tmp/reduce.mpileup'
    
    # Delete previous run files if they exist
    if(os.path.exists(temp_mpileup)):
        os.remove(temp_mpileup)

    # S3 SELECT query to get the rows where the second column is in the selected range
    expression = "SELECT * FROM s3object s WHERE cast(s._2 as int) BETWEEN %s AND %s" % (range['start'], range['end'])
    input_serialization = {'CSV': {'RecordDelimiter': '\n', 'FieldDelimiter': '\t'}, 'CompressionType': 'NONE'}

    # Execute S3 SELECT
    stat.timer_start(expression)
    for k in keys:
        subStat.timer_start(k)
        try:
            resp = s3.select_object_content(
                Bucket=pipeline_params.storage_bucket,
                Key=k,
                ExpressionType='SQL',
                Expression=expression,
                InputSerialization = input_serialization,
                OutputSerialization = {'CSV': {"FieldDelimiter" : "\t"}}
            )
        except:
            raise ValueError("ERROR IN KEY: " + k)

        data = ""
        for event in resp['Payload']:
            if 'Records' in event:
                records = event['Records']['Payload'].decode("UTF-8")
                data = data + records

        with open(temp_mpileup, 'a') as f:
            f.write(data)
        del data
        subStat.timer_stop(k)
    stat.timer_stop(expression)
    stat.store_dictio(subStat.get_stats(), "subprocesses", expression)

    # Execute the script to merge and reduce
    sinple_out = sp.check_output(['bash', '/function/bin/mpileup_merge_reducev3_nosinple.sh', temp_mpileup, '/function/bin/', "75%"])
    sinple_out = sinple_out.decode('UTF-8')

    # Final file
    sinple_name=temp_mpileup+'_merged.mpileup'

    # write output to /tmp
    with open(sinple_name, 'w') as f:
        f.write(sinple_out)
    
    #Upload part
    part = s3.upload_part(
        Body = sinple_out,
        Bucket = pipeline_params.storage_bucket,
        Key = mpu_key,
        UploadId = mpu_id,
        PartNumber = n_part
    )
    
    mainStat.timer_stop("call_reduceFunction")
    mainStat.store_dictio(stat.get_stats(), "subprocesses", "call_reduceFunction")

    return {"PartNumber" : n_part, "ETag" : part["ETag"], "mpu_id": mpu_id}, mainStat.get_stats()