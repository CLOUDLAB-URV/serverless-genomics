1st, install python 3.9

Set up dependencies:

```
python3.9 -m venv venv && source venv/bin/activate && pip install -r requirements.txt
```

Runtime docker images have to be built and uploaded to ECR:

AWS Lambda:
```console
lithops runtime build -f Dockerfile-aiobotocore-lambda -b aws_lambda aiobotocore-vanilla-pyv39:v0
```

AWS Batch:
```console
lithops runtime build -b aws_batch -f Dockerfile-aiobotocore-batch batch_aiobotocore:v0
```

Then, the lithops wrapper can be executed like this (--memory is optional):

```console
python LaunchJob.py --memory 1770 bucketname inputs/SRR6052133_2_split_1.se.map -- ./MergeOrPair.native merge -t 2
```

.lithops config:

```
aws_s3:
    storage_bucket: gil-lithops-meta-cloudbutton-hutton
    region_name : us-east-1

aws_batch:
    execution_role: arn:aws:iam::423521657208:role/ecsTaskExecutionRole
    instance_role: arn:aws:iam::423521657208:instance-profile/ecsInstanceRole
    region_name : us-east-1
    env_type: SPOT
    env_max_cpus: 2
    container_vcpus: 2
    #assign_public_ip: no
    subnets:
        - subnet-0c2b949cd2d35acbe
    security_groups:
        - sg-0cac347c598f40ea5
    runtime_memory: 3072
    runtime_timeout: 3600
```
