#!/bin/bash
#
# NOTE: you must have aws-cli installed and configured to use this
# ALSO NOTE: This will use the default VPC/ subnets AWS creates for you in your specified region
#
REGION="us-east-1"
PEM="tyler_key.pem"
KEY="tyler_key"
SECURITY_GROUP="tyler_group"
FORCE_GROUP_PERMISSIONS=0
LOG="log.txt"
EBS_SIZE=

echo -e "\nStarting up new EC2 instance...\n"

# set region
aws configure set default.region $REGION

#generate key pair if provided PEM file doesn't exist
if [ -f $PEM ]
then
  echo -e "$PEM exists, skipping key-pair creation\n" | tee -a $LOG
else
  echo -e "No PEM file provided, generating key-pair with ID $KEY\n" | tee -a $LOG
  res=`aws ec2 create-key-pair --key-name $KEY`
  echo $res | tee -a $LOG
  key=`echo $res | python3 -c "import sys, json; print(json.load(sys.stdin)['KeyMaterial'])"`
  PEM=$key".pem"
  echo $key > $PEM
fi
PEM=$(realpath "$PEM")


#generate security group
sgroup_json=`aws ec2 describe-security-groups --filters Name=group-name,Values=$SECURITY_GROUP`
GROUP_ID="NAN"
group_exists=0
if echo $sgroup_json | grep -q $SECURITY_GROUP
then
  group_exists=1
  echo -e "Security group $SECURITY_GROUP exists... Skipping security group configuration.\n"  | tee -a $LOG
  #echo $sgroup_json
  GROUP_ID=`echo $sgroup_json | python3 -c "import sys, json; print(json.load(sys.stdin)['SecurityGroups'][0]['GroupId'])"`

else
  ret=`aws ec2 create-security-group --group-name $SECURITY_GROUP --description "Security group for lithops variant_caller"`
  echo $ret | tee -a  $LOG
  GROUP_ID=`echo $ret | python3 -c "import sys, json; print(json.load(sys.stdin)['GroupId'])"`
fi
echo -e "Fetched Security Group ID: $GROUP_ID\n" | tee -a $LOG

# make sure security ID is valid
if [ "$GROUP_ID" == "NAN" ] || [ "$GROUP_ID" == "" ]
then
  echo -e "Invalid group id\n"  | tee -a $LOG
  exit 1
fi

# set group security rules
to_set=0
if [ $group_exists -eq 1 ]
then
  if [ $FORCE_GROUP_PERMISSIONS -eq 1 ]
  then
    echo "Forcing configuration of existing group $SECURITY_GROUP (FORCE_GROUP_PERMISSIONS==1)..."  | tee -a $LOG
    to_set=1
  fi
else
  to_set=1
fi
if [ $to_set -eq 1 ]
then
  echo "Setting security group permissions for $GROUP_ID"  | tee -a $LOG

  # define inbound traffic rules
  # note we allow some thing here that aren't necessary for the client VM, this
  #   is so we can re-use the same security group for the redis server VM. If
  #   Not desired, comment out the redis and infiniband lines below
  #ssh inbound
  aws ec2 authorize-security-group-ingress --region $REGION --group-id $GROUP_ID --protocol tcp --port 22 --cidr '0.0.0.0/0' | tee -a  $LOG
  # redis inbound
  aws ec2 authorize-security-group-ingress --region $REGION --group-id $GROUP_ID --protocol tcp --port 6379 --cidr '0.0.0.0/0' | tee -a  $LOG
  # infiniband inbound
  aws ec2 authorize-security-group-ingress --region $REGION --group-id $GROUP_ID --protocol tcp --port 11222 --cidr '0.0.0.0/0' | tee -a  $LOG
  # docker daemon
  aws ec2 authorize-security-group-ingress --region $REGION --group-id $GROUP_ID --protocol tcp --port 2375 --cidr '0.0.0.0/0' | tee -a  $LOG
  # docker daemon encrypted
  aws ec2 authorize-security-group-ingress --region $REGION --group-id $GROUP_ID --protocol tcp --port 2376 --cidr '0.0.0.0/0' | tee -a  $LOG

  # allow all outbound
  aws ec2 authorize-security-group-egress --region $REGION --group-id $GROUP_ID --protocol all | tee -a  $LOG
fi

# get latest ubuntu 20.04 AMI
images=`aws ec2 describe-images --filters "Name=name,Values=ubuntu*20.04-arm64-server*" --query "sort_by(Images, &CreationDate)[-1:]"`
IMAGE_NAME=`echo $images | python3 -c "import sys, json; print(json.load(sys.stdin)[0]['Name'])"`
AMI=`echo $images | python3 -c "import sys, json; print(json.load(sys.stdin)[0]['ImageId'])"`
echo -e "Using image $IMAGE_NAME with AMI $AMI\n" | tee -a  $LOG


#aws ec2 run-instances --image-id ami-xxxxxxxx --enable-api-termination --security-group-ids xxxxxxxxx --instance-type t2.micro --key-name example-key
INSTANCE=tt

# Output report
echo -e "
#######################################################################
Done! Created running instance with the following parameters:

REGION=$REGION
PEM=$PEM
KEY=$KEY
SECURITY_GROUP=$SECURITY_GROUP
GROUP_ID=$GROUP_ID
FORCE_GROUP_PERMISSIONS=$FORCE_GROUP_PERMISSIONS
IMAGE_NAME=$IMAGE_NAME
AMI=$AMI
INSTANCE_NAME=$INSTANCE_NAME
INSTANCE=$INSTANCE

All outputs are saved in $LOG
"

# output instructions to user to shut down instance
echo -e "
#######################################################################
Note that your AWS account will be charged for running instances.

To check for running instances, use:
aws ec2 describe-instances

To check on this specific instance, use:
aws ec2 describe-instances --instance-ids $INSTANCE

To **stop** this instance, use:
aws ec2 stop-instances --instance-ids $INSTANCE

[[Note that AWS does not charge for **stopped** instances, however you may still be charged for EBS volume usage.]]

To **terminate** this instance, use:
aws ec2 terminate-instances --instance-ids $INSTANCE
" | tee -a  $LOG
