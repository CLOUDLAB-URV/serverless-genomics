import sys
import os
import time
import re

fq_set_n=5
keys=('/tmp/ERR9729866_fq12-fa0_map.index.txt', 
      '/tmp/ERR9729866_fq12-fa1_map.index.txt',
      '/tmp/ERR9729866_fq12-fa2_map.index.txt',
      '/tmp/ERR9729866_fq12-fa3_map.index.txt',
      '/tmp/ERR9729866_fq12-fa4_map.index.txt',
      '/tmp/ERR9729866_fq6-fa0_map.index.txt',
      '/tmp/ERR9729866_fq6-fa1_map.index.txt',
      '/tmp/ERR9729866_fq6-fa2_map.index.txt',
      '/tmp/ERR9729866_fq6-fa3_map.index.txt',
      '/tmp/ERR9729866_fq6-fa4_map.index.txt',
      '/tmp/ERR9729866_fq3-fa1_map.index.txt',
      '/tmp/ERR9729866_fq3-fa2_map.index.txt',)
fq_set_dict={}

for key in keys:
    fq_set=re.sub(r'/tmp/(\S*_fq\d*)-fa\d*_map.index.txt', r'\1', key)
    file_name=re.sub(r'/tmp/','',key)

    #print("fq set: "+str(fq_set))
    #print("file name: "+file_name)
    
    # download file

    # add key to fq_set dictionary and count number of occurrences of key
    if fq_set in fq_set_dict:
        fq_set_dict[fq_set] = int(fq_set_dict[fq_set]) + 1
    else:
        fq_set_dict[fq_set] = 1

# retrieve ready sets
ready_sets=[]
for key, value in fq_set_dict.items():
    print(key+ "\t" + str(value))
    if value == 5:
        ready_sets.append(key)
        

# complete sets
print("\nready sets")
for set in ready_sets:
    print(set)
    # delete key of complete sets
    del fq_set_dict[set]

print("\nremaining sets")
for key, value in fq_set_dict.items():
    print(key+ "\t" + str(value))
