import aws_ec2_control as ec2
import time

EC2_IP = '54.146.89.181'
EC2_ID = 'i-0cee52f66655d990b'
REGION_NAME = 'us-east-1'

#actions:
# 1. start server, clean redis and stop
# 2. clean redis and stop
# 3. start server and clean redis
# 4. clean redis

action=4

if action==1:
    ec2.action_on_instance(region_name=REGION_NAME,instance_id=EC2_ID,action='start')
    time.sleep(25)
    ec2.flush_redis_database(ec2_IP_address=EC2_IP,redis_password='lucio-redis-server')
    #time.sleep(10)
    #ec2.start_redis_reducer(path_to_secret_key='~/bioss/cloudbutton/variant_caller/v3/aws-redis-damien.pem',ec2_IP_address=EC2_IP,iteration_number=10)
    time.sleep(5)
    ec2.action_on_instance(region_name=REGION_NAME,instance_id=EC2_ID,action='stop')
elif action==2:
    ec2.flush_redis_database(ec2_IP_address=EC2_IP,redis_password='lucio-redis-server')
    time.sleep(5)
    ec2.action_on_instance(region_name=REGION_NAME,instance_id=EC2_ID,action='stop')
elif action==3:
    ec2.action_on_instance(region_name=REGION_NAME,instance_id=EC2_ID,action='start')
    time.sleep(25)
    ec2.flush_redis_database(ec2_IP_address=EC2_IP,redis_password='lucio-redis-server')
elif action==4:
    ec2.flush_redis_database(ec2_IP_address=EC2_IP,redis_password='lucio-redis-server')
else:
    print("no ec2 action specified")
