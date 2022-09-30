import lithops
import os
import subprocess as sp
from lithops import Storage

file_list = []

def download_files_to_local(bucket_name,redis):
    keys = redis.list_keys(bucket_name)
    #new_keys = {x.replace('/tmp','').replace('/tmp','').replace('.gem','.sam') for x in keys}
    new_keys = {x.replace('/tmp/','') for x in keys}
    for files, new_files in zip(keys,new_keys):
        obj_stream = redis.get_object(bucket_name, files)
        file = open(new_files, 'wb')
        file.write(obj_stream)
        file.close()
    return (keys, new_keys)

def group_by_fasta(file_list):
    splits = []
    for i in file_list:
        splits.append(i[:14])
    splits = list(dict.fromkeys(splits))
    return(splits)

def execute_binary_reduce(splits,redis):
    for split in splits:
        print('./binary_reducer merge_gem_alignment_metrics.sh 1 '+split+'* >  '+split+'.intermediate.txt')
        cmd = f'./binary_reducer.sh /home/ubuntu/merge_gem_alignment_metrics.sh 1 {split}* > {split}.intermediate.txt'
        sp.run(cmd, shell=True, check=True, universal_newlines=True)
        print('./filtered_merged_index.sh '+split+'.intermediate.txt '+split)
        cmd2 = f'./filter_merged_index.sh {split}.intermediate.txt {split}'
        sp.run(cmd2, shell=True, check=True, universal_newlines=True)
        with open(split+'.txt','r') as f:
            redis.put_object('',split+'.txt',f.read())
    return('Binary reduce executed!')

def delete_files_from_local(file_list,splits):
    for file in file_list:
        os.remove(file)
    for split in splits:
        os.remove(split+'.intermediate.txt')
        os.remove(split+'.txt')
    return('Temporary local files removed!')

def delete_files_from_redis(keys_list,redis):
    redis.delete_objects('',keys_list)
    return('Files deleted from redis!')


if __name__ == "__main__":
   redis = Storage(backend='redis')
   keys_list, file_list = download_files_to_local('',redis)
   splits = group_by_fasta(file_list)
   print(splits)
   print(execute_binary_reduce(splits,redis))
   print(delete_files_from_local(file_list,splits))
   print(delete_files_from_redis(keys_list,redis))