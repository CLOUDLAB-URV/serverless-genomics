

# 1. ON LAPTOP
# update credentials
scp -i "aws-redis-damien.pem" /home/lumar/.lithops/aws_credentials ubuntu@ec2-34-236-249-86.compute-1.amazonaws.com:/home/ubuntu/lucio
scp -i "aws-redis-damien.pem" /home/lumar/.lithops/lithops_config ubuntu@ec2-34-236-249-86.compute-1.amazonaws.com:/home/ubuntu/lucio

# connect to ec2
ssh -i aws-redis-damien.pem ubuntu@ec2-34-236-249-86.compute-1.amazonaws.com

# update variant caller, varcall_args and map reduce
scp -i "aws-redis-damien.pem" varcall_lithops_demo_v7.py ubuntu@ec2-34-236-249-86.compute-1.amazonaws.com:/home/ubuntu/lucio/serverless_genomics/variant_caller
scp -i "aws-redis-damien.pem" varcall_args.tsv ubuntu@ec2-34-236-249-86.compute-1.amazonaws.com:/home/ubuntu/lucio/serverless_genomics/variant_caller
scp -i "aws-redis-damien.pem" map_reduce.py ubuntu@ec2-34-236-249-86.compute-1.amazonaws.com:/home/ubuntu/lucio/serverless_genomics/variant_caller

# retrieve varcall_args
scp -i "aws-redis-damien.pem"  ubuntu@ec2-34-236-249-86.compute-1.amazonaws.com:/home/ubuntu/lucio/serverless_genomics/variant_caller/varcall_args*  varcall_out/
# retrieve varcall logs
scp -i "aws-redis-damien.pem"  ubuntu@ec2-34-236-249-86.compute-1.amazonaws.com:/home/ubuntu/lucio/serverless_genomics/variant_caller/varcall_out/varcall_89_1_*rep1.log  varcall_out/
scp -i "aws-redis-damien.pem"  ubuntu@ec2-34-236-249-86.compute-1.amazonaws.com:/home/ubuntu/lucio/serverless_genomics/variant_caller/varcall_out/varcall_92_2_*  varcall_out/
scp -i "aws-redis-damien.pem"  ubuntu@ec2-34-236-249-86.compute-1.amazonaws.com:/home/ubuntu/lucio/serverless_genomics/variant_caller/varcall_out/*  varcall_out/docker


##################
# 2. ON EC2
# run docker image using convenience script 
cd lucio
./serverless_genomics/variant_caller/varcall_client/varcall_docker_run.sh ./serverless_genomics/ aws_credentials lithops_config tkchafin/varcall_client:0.2
cd variant_caller
# install unbuffer in dockerfile
apt-get -y install expect
#check active run details
grep -v "^#" varcall_args.tsv 

bash run_variant_caller_v2.sh varcall_args.tsv 2 1 True


rm varcall_lithops_demo_v7.py
rm varcall_args.tsv
rm map_reduce.py

#pull code changes to ec2 with dockerfile - over-writing local files
git fetch --all
git reset --hard origin/main


# Parse information from latest log
head -n 30 $(ls -tr *rep1.log | tail -n 1)
tail -n 60 $(ls -tr *rep1.log | tail -n 1)
more  $(ls -tr *rep1.log | tail -n 1)
wc -l $(ls -tr *rep1.log | tail -n 1)
head -n 547350 $(ls -tr *rep1.log | tail -n 1) | tail -n 100

# get iterdata number
grep "fasta x fastq chunks"  $(ls -tr *rep1.log | tail -n 1) 
# gem mapper
awk 'NF==11 && match($11, /filt_wline_no.map$/) {if (!a[$0]++) {print $0}}' $(ls -tr *rep1.log | tail -n 1) | sort --uniq | wc -l
# gempileup
grep "_filt_wline_no_corrected.map.mpileup$"  $(ls -tr *rep1.log | tail -n 1) | sort --uniq | wc -l
# conversion to csv
grep "filt_wline_no_corrected.map.mpileup.csv" $(ls -tr *rep1.log | tail -n 1) | sort --uniq| wc -l
# initiation / completion of map / reduce functios
grep -n "Getting results from" $(ls -tr *rep1.log | tail -n 1)
# look at map function lithops log
grep "ExecutorID 5e6771-0" $(ls -tr *rep1.log | tail -n 1)
# map function completion
grep "Pending: 0 - Running: 0" $(ls -tr *rep1.log | tail -n 1)

# reduce phase
grep "reduce.mpileup_merged.mpileup$" $(ls -tr *rep1.log | tail -n 1) | sort --uniq | wc -l
grep -n "count indexes: execution_time" $(ls -tr *rep1.log | tail -n 1)
grep -n "s3 select: execution_time" $(ls -tr *rep1.log | tail -n 1)

grep "ExecutorID 5e6771" $(ls -tr *rep1.log | tail -n 1) 
grep "ExecutorID 5e6771" $(ls -tr *rep1.log | tail -n 1) | wc -l
grep "ExecutorID 5e6771" $(ls -tr *rep1.log | tail -n 1) | head -n 20
grep "JobID M000 - Calls" $(ls -tr *rep1.log | tail -n 1) | wc -l
grep "ExecutorID 5e6771 | JobID M000 - Calls"  $(ls -tr *rep1.log | tail -n 1) | wc -l