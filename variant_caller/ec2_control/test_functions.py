import aws_ec2_control as ec2

EC2_IP = '54.146.89.181'
EC2_ID = 'i-0cee52f66655d990b'
REGION_NAME = 'us-east-1'


#ec2.action_on_instance(region_name=REGION_NAME,instance_id=EC2_ID,action='start')
#
ec2.action_on_instance(region_name=REGION_NAME,instance_id=EC2_ID,action='stop')
#
#ec2.flush_redis_database(ec2_IP_address=EC2_IP,redis_password='lucio-redis-server')
#
#ec2.start_redis_reducer(path_to_secret_key='~/bioss/cloudbutton/variant_caller/v3/aws-redis-damien.pem',ec2_IP_address=EC2_IP,iteration_number=10)
