#
# (C) Copyright Cloudlab URV 2021
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

import boto3
import subprocess as sp



def get_EC2_con(my_region):
    get_EC2_con_re=boto3.resource('ec2',region_name=my_region)
    return get_EC2_con_re

def list_instances_in_my_region(ec2_con_re):
    for each in ec2_con_re.instances.all():
        print (each.id)

def get_instant_state(ec2_con_re,in_id):
    for each in ec2_con_re.instances.filter(Filters=[{'Name':'instance-id','Values':[in_id]}]):
        pr_st=each.state['Name']
    return pr_st

def start_instance(ec2_con_re, in_id):
    pr_st = get_instant_state(ec2_con_re,in_id)
    if pr_st == "running":
        print("Instance is already running")
    else:
        for each in ec2_con_re.instances.filter(Filters=[{'Name':'instance-id','Values':[in_id]}]):
            each.start()
            print("Please wait until the EC2 starts...")
            each.wait_until_running()
            print("Running!")
    return

def stop_instance(ec2_con_re,in_id):
    pr_st = get_instant_state(ec2_con_re,in_id)
    if pr_st == "stopped":
        print("Instance is already stopped")
    else:
        for each in ec2_con_re.instances.filter(Filters=[{'Name':'instance-id','Values':[in_id]}]):
            each.stop()
            print("Please wait until the EC2 stops...")
            each.wait_until_running()
            print("Stopped!")
    return

def action_on_instance(region_name,instance_id,action):
    ec2_con_re=get_EC2_con(region_name)
    while True:
        if action not in ["start","stop"]:
            print("Please enter an accepted action (start or stop)")
        else:
            break
    if action == "start":
        start_instance(ec2_con_re,instance_id)
    else:
        stop_instance(ec2_con_re,instance_id)

def flush_redis_database(path_to_secret_key,ec2_IP_address,redis_password):
    print("Flushing the redis database")
    cmd = f'redis-cli -h {ec2_IP_address} -p 6379 -a {redis_password} FLUSHDB'
    sp.run(cmd, shell=True, check=True, universal_newlines=True)
    # kill redis_server.py it is hasn't stopped
    cmd = f"ssh -o StrictHostKeyChecking=no -i {path_to_secret_key} ubuntu@{ec2_IP_address} 'pkill -9 -f redis_server.py'"
    sp.Popen(cmd, shell=True, universal_newlines=True)

def start_redis_reducer(path_to_secret_key,ec2_IP_address, fastq_set_n, run_id, iterdata_map_n):
    cmd = f"ssh -o StrictHostKeyChecking=no -i {path_to_secret_key} ubuntu@{ec2_IP_address} 'python3 redis_server.py {fastq_set_n} {run_id} {iterdata_map_n} &'"
    #sp.run(cmd, shell=True, check=True, universal_newlines=True)
    sp.Popen(cmd, shell=True, universal_newlines=True)
