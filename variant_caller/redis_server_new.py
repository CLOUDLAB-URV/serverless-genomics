# ssh -i "aws-redis-damien.pem" ubuntu@ec2-54-146-89-181.compute-1.amazonaws.com
# scp -i "aws-redis-damien.pem" ubuntu@ec2-54-146-89-181.compute-1.amazonaws.com:/home/ubuntu/redis_server.py .
# scp -i "aws-redis-damien.pem" ubuntu@ec2-54-146-89-181.compute-1.amazonaws.com:/home/ubuntu/redis_lm.log .
# scp -i "aws-redis-damien.pem" redis_server.py ubuntu@ec2-54-146-89-181.compute-1.amazonaws.com:/home/ubuntu/ 
import sys
import os
import re
import time
import subprocess as sp
import redis

def process_redis_keys(r,fq_set_dict,key_dict,key_count):
    
    ready_sets=[]
    files_to_delete=[]
    
    keys = r.keys("*")
    if len(keys) >= 1 and len(keys[0]) > 10:
        
        log.write("\nlisting redis keys\n")
        log.write("keys type: "+str(type(keys))+"\n")
        
        log.write("\n#NEW "+ run_id +" ITERATION\nkeys: "+str(keys)+"\n")
        for key in keys:
            key = key.decode('utf-8')
            print(run_id,key)
            print(key.find(str(run_id)))
            
            #If the key read from the variant caller is a duplicate, continue to the next element
            if key in key_dict:
                    log.write("duplicate key: "+str(key) +"\n")
                    #r.delete(key)
                    continue
            #Else, process the new key
            else:
                if key.find(str(run_id)) != -1:
                    print("key found: "+key)
                    # obtain fq set and file name from key
                    fq_set=re.sub(r'/tmp/(\S*_fq\d*)-fa\d*_map.index.txt', r'\1', key)
                    file_name=re.sub(r'/tmp/','',key)
                    remove_run_id=run_id+"__"
                    fq_set=re.sub(remove_run_id,'',fq_set)
                    file_name=re.sub(remove_run_id,'',file_name)

                    key_count+=1
                    key_dict[key] = key
                    log.write("key "+str(key_count)+": "+str(key) +"\n")
                
                    # download file
                    obj_stream = r.get(key)
                    file = open(file_name, 'wb')
                    file.write(obj_stream)
                    file.close()

                    # add key to fq_set dictionary and count number of occurrences of key
                    if fq_set in fq_set_dict:
                        fq_set_dict[fq_set] = int(fq_set_dict[fq_set]) + 1
                    else:
                            fq_set_dict[fq_set] = 1
                else:
                    log.write("sporious key: "+str(key) +"\n")
                
                #print("fq set: "+str(fq_set))
                #print("file name: "+file_name)
                
            



        # identify complete sets in dictionary 
        keys_in_fq_set_dict=0
        log.write("fastq set dictionary status:\n")
        for key, value in fq_set_dict.items():
            #log.write('key: ' + str(key) + ' fq_set_n: '+ fq_set_n + "\t" + ' value:'+ str(value)+"\n")
            log.write('val and fq: '+ str(value) + ' ' + str(fq_set_n))
            if value == fq_set_n:
                log.write("found complete set: "+key+"\n")
                ready_sets.append(key)
            keys_in_fq_set_dict= keys_in_fq_set_dict + value
        log.write("number of keys in fq set dictionary: " + str(keys_in_fq_set_dict)+"\n")
        log.write("number of keys in keys dictionary:   " + str(len(key_dict))+"\n")
        for set in ready_sets:
            for i in range(0, fq_set_n-1):
                file_to_delete=str(set)+"-fa"+str(i)+"_map.index.txt"
                if file_to_delete in key_dict:
                    files_to_delete.append(file_to_delete)
                else:
                    print(file_to_delete+" not found in keys dictionary - incomplete set")

        # delete keys from redis
        #r.delete(*keys) 
        
        

    return (ready_sets, files_to_delete, key_count)


def execute_binary_reduce(set,r):
    log.write('./binary_reducer merge_gem_alignment_metrics.sh 1 '+set+'* >  '+set+'.intermediate.txt'+"\n")
    cmd = f'./binary_reducer.sh ./merge_gem_alignment_metrics.sh 1 {set}* > {set}.intermediate.txt'
    print(cmd)
    sp.run(cmd, shell=True, check=True, universal_newlines=True)
    log.write('./filtered_merged_index.sh '+set+'.intermediate.txt '+set+"\n")
    cmd2 = f'./filter_merged_index.sh {set}.intermediate.txt {set}'
    print(cmd2)
    sp.run(cmd2, shell=True, check=True, universal_newlines=True)
    with open(set+'.txt','r') as f:
        log.write('put r:' + str(set))
        r.set(set+'.txt',f.read())
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
    
    fq_set_n = int(sys.argv[1])
    run_id = sys.argv[2]
    iterdata_n = int(sys.argv[3])

    log.write("\n\nRUN "+run_id+" START: redis server index processing - "+time.strftime("%Y%m%d_%H:%M:%S")+"\n")
    log.write("iterdata number: "+ str(iterdata_n)+"\n")
    log.write("fastq set no.: "+ str(fq_set_n)+"\n")

    r = redis.Redis(host='localhost', port=6379,password='lucio-redis-server',socket_connect_timeout=1)

    # dictionary for all incoming keys
    key_dict={}
    fq_set_dict={}
    key_count=0
    while (True):
        log.write("\n##### "+run_id+" MAIN ITERATION\n")
        ready_sets, files_to_delete, key_count = process_redis_keys( r, fq_set_dict, key_dict, key_count)
        log.write("type of ready_sets: "+str(type(ready_sets)))
        if len(ready_sets) >= 1:
            log.write("processing ready sets: "+str(ready_sets)+"\n")
            log.write("files to be deleted: "+str(files_to_delete)+"\n")
            log.write("key count: "+str(key_count)+"\n")
            for set in ready_sets:
                execute_binary_reduce(set,r)
                # delete key of complete sets
                del fq_set_dict[set]
        #delete_files_from_local(files_to_delete , ready_sets)
        log.flush()    
        if int(iterdata_n) == key_count:
            break
        time.sleep(1)
