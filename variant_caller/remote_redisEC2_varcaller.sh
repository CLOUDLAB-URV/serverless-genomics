#!/bin/bash

execute=false
if [ "$execute" = true ]; then
    ssh -i "aws-redis-damien.pem" ubuntu@ec2-54-146-89-181.compute-1.amazonaws.com
fi
# kill any redis_server.py scripts still running
ssh -i "aws-redis-damien.pem" ubuntu@ec2-54-146-89-181.compute-1.amazonaws.com 'pkill -9 -f redis_server.py'
# update redis.server.py (local to VM)
scp -i "aws-redis-damien.pem" redis_server.py ubuntu@ec2-54-146-89-181.compute-1.amazonaws.com:/home/ubuntu/
# find txt filenames and append to log
ssh -i "aws-redis-damien.pem" ubuntu@ec2-54-146-89-181.compute-1.amazonaws.com 'find . -maxdepth 1 -name "*txt" | sort >> redis_lm.log &'
# retrieve log (VM to local)
scp -i "aws-redis-damien.pem" ubuntu@ec2-54-146-89-181.compute-1.amazonaws.com:/home/ubuntu/redis_lm.log  "redis_lm_"$(date +%Y%m%d_%H%M)".log"
# delete text filenames
ssh -i "aws-redis-damien.pem" ubuntu@ec2-54-146-89-181.compute-1.amazonaws.com 'rm *txt *log &'
# clear redis 
python ec2_control/redis_maintenance.py 

scp -i "aws-redis-damien.pem" ubuntu@ec2-54-146-89-181.compute-1.amazonaws.com:/home/ubuntu/filter_merged_index.sh .

scp -i "aws-redis-damien.pem" ubuntu@ec2-54-146-89-181.compute-1.amazonaws.com:/home/ubuntu/redis_server.py  "redis_server_new.py"