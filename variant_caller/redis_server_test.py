# ssh -i "aws-redis-damien.pem" ubuntu@ec2-54-146-89-181.compute-1.amazonaws.com
# scp -i "aws-redis-damien.pem" ubuntu@ec2-54-146-89-181.compute-1.amazonaws.com:/home/ubuntu/redis_server.py .
# scp -i "aws-redis-damien.pem" ubuntu@ec2-54-146-89-181.compute-1.amazonaws.com:/home/ubuntu/redis_lm.log .
# scp -i "aws-redis-damien.pem" redis_server.py ubuntu@ec2-54-146-89-181.compute-1.amazonaws.com:/home/ubuntu/ 
import lithops
import sys
import os
import re
import time
import subprocess as sp
from lithops import Storage


def process_redis_keys(bucket_name,redis,fq_set_dict,key_count):
    
    ready_sets=[]
    files_to_delete=[]
    
    keys = redis.list_keys(bucket_name)
    if len(keys) >= 1:
        
        #log.write("\nlisting redis keys\n")
        #log.write("keys type: "+str(type(keys))+"\n")
        log.write("\nkeys: "+str(keys)+"\n")

        for key in keys:
            key_count+=1
            fq_set=re.sub(r'/tmp/(\S*_fq\d*)-fa\d*_map.index.txt', r'\1', key)
            file_name=re.sub(r'/tmp/','',key)
            log.write("key "+str(key_count)+": "+str(key) +"\n")
            #print("fq set: "+str(fq_set))
            #print("file name: "+file_name)
            
            # download file
            obj_stream = redis.get_object(bucket_name, key)
            file = open(file_name, 'wb')
            file.write(obj_stream)
            file.close()

            # add key to fq_set dictionary and count number of occurrences of key
            if fq_set in fq_set_dict:
                fq_set_dict[fq_set] = int(fq_set_dict[fq_set]) + 1
            else:
                fq_set_dict[fq_set] = 1

        # identify complete sets in dictionary 
        ready_sets=[]
        for key, value in fq_set_dict.items():
            print(key+ "\t" + str(value))
            if value == fq_set_n:
                log.write("found complete set: "+key+"\n")
                ready_sets.append(key)
        files_to_delete=[]
        for set in ready_sets:
            for i in range(0, fq_set_n-1):
                file_to_delete=set+"-fa"+i+"_map.index.txt"
                files_to_delete.append(file_to_delete)

        # delete keys from redis
        redis.delete_objects('',keys) 
        

    return (ready_sets, files_to_delete, key_count)


def execute_binary_reduce(set,redis):
    log.write('./binary_reducer merge_gem_alignment_metrics.sh 4 '+set+'* >  '+set+'.intermediate.txt'+"\n")
    cmd = f'./binary_reducer.sh /home/ubuntu/merge_gem_alignment_metrics.sh 4 {set}* > {set}.intermediate.txt'
    sp.run(cmd, shell=True, check=True, universal_newlines=True)
    log.write('./filtered_merged_index.sh '+set+'.intermediate.txt '+set+"\n")
    cmd2 = f'./filter_merged_index.sh {set}.intermediate.txt {set}'
    sp.run(cmd2, shell=True, check=True, universal_newlines=True)
    with open(set+'.txt','r') as f:
        redis.put_object('',set+'.txt',f.read())
    return('Binary reduce executed!')

def delete_files_from_local(file_list,sets):
    for file in file_list:
        os.remove(file)

    for set in sets:
        os.remove(set+'.intermediate.txt')
        os.remove(set+'.txt')
    return('Temporary local files removed!')

if __name__ == "__main__":
    log=open("redis_lm.log", "w") 
    log.write("starting redis server index processing\n")
    fq_set_n = int(sys.argv[1])
    log.write("fastq set no.: "+ str(fq_set_n)+"\n")
    fq_set_dict={}
    key_count=0
    while (True):
        key_count+=1
        if key_count<50:
            log.write("loop\n")
        #time.sleep(0.3)
    