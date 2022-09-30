# Running Docker in Amazon EC2

Here, we provide instructions to run the variant_caller client using a pre-built Docker image running in an AWS EC2 VM. Docker is a platform that allows you to build and deploy packaged software environments ("containers"). Combining this with the use of AWS Virtual Machines (VMs), you can get up-and-running quickly, without the need to leverage local resources.

## 1. Virtual Machine configuration

### 1.1 VM setup on AWS

Configuring the virtual machine instance follow the same process used for configuring the redis server, with only a small exception with respect to the storage configuration (see below). Note that the below steps have only been tested on the Ubuntu 20.04 instance type, which is free-tier eligible, thus we recommend this as the OS selection.

To set up your AWS account for running the lithops variant caller in a VM, see [this guide](https://github.com/Damian-MG/CloudButton-Redis-Installation/blob/main/AWS_account_configuration.md).

Once your account is setup, reference [section 1 (Virtual Machine Configuration) from the redis-server setup documentation](https://github.com/Damian-MG/CloudButton-Redis-Installation).

### 1.2 Storage configuration

When selecting your disk storage size, the total size will at minimum need to be large enough to host the docker image ([tkchafin/varcall_client](https://hub.docker.com/repository/docker/tkchafin/varcall_client)). Uncompressed, this is >3GB. Additionally, sufficient space will be needed to temporarily store the input FASTQ file in /tmp during the pre-processing phase, unless indexes are already present in the S3 bucket (see [main pipeline documentation]()).

### 1.3 Connecting to the VM via SSH

If configured correctly to allow incoming SSH traffic (see [section 1.4: Network Configuration from the redis-server documentation](https://github.com/Damian-MG/CloudButton-Redis-Installation)), you should be able to now SSH directly to the running instance using the private key file specified in the VM configuration:
`ssh -i private_key.pem ubuntu@<VirtualMachine_DNS>`

Note that the username is `ubuntu` and the <VirtualMachine_DNS> may be found in the instance summary, under "Public IPv4 DNS" (it will look like `ec2-##-###-##-##.compute-1.amazonaws.com`).

After running this, you should notice your terminal prompt has changed to:
```
ubuntu@ip-###-##-##-##:~$
```

## 2. Docker configuration in the VM

Once you have verified that you can successfully access the VM instance from your local machine, you are ready to configure the docker. Note that here there are several steps which will need to be completed either from the your local computer, and others from the VM, accessed using SSH (see Section 1.3 above).

For convenience, we recommend opening two terminals, one running locally, and another connected to the VM via SSH.  

### 2.1 Installing docker in the VM

This step occurs on the **VM**. Here, you will install docker. Docker installation follows the [Docker documentation](https://docs.docker.com/engine/install/ubuntu/), with commands excerpted below for convenience (with some slight modifications):

```
# ssh to the running VM instance
ssh -i private_key.pem ubuntu@<VirtualMachine_DNS>

# remove previous versions (just in case)
sudo apt-get remove docker docker-engine docker.io containerd runc

# set up the repository
sudo apt-get update && sudo apt-get install -y ca-certificates curl gnupg lsb-release

# add Docker's official GPG key
sudo mkdir -p /etc/apt/keyrings
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /etc/apt/keyrings/docker.gpg

# set up repo
echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.gpg] https://download.docker.com/linux/ubuntu \
  $(lsb_release -cs) stable" | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null

# install docker engine
sudo apt-get update && sudo apt-get install -y docker-ce docker-ce-cli containerd.io docker-compose-plugin

# add user to docker group
sudo usermod -a -G docker ubuntu

```

Now, you can pull the pre-built docker image from our public DockerHub repository:

```
# docker pull
docker pull tkchafin/varcall_client:0.2

```

If you want to make modifications to the docker image, or build it yourself, see Section 3: Building the container from scratch (below).

### 2.2 Install variant_caller code on the VM

Next, you will need a copy of the variant_caller repository on the **VM**. From your **VM** SSH terminal, first install git and then clone the repo:
```
# install git
sudo apt-get install -y git

# clone repository
git clone https://gitlab1.bioss.ac.uk/lmarcello/serverless_genomics.git

# if using a development branch, checkout that branch, e.g.,
# git checkout some_dev_branch

```

### 2.3 Upload user configuration files to VM

Now, from your **local** terminal, you need to upload your [lithops configuration](https://lithops-cloud.github.io/docs/source/configuration.html) and [AWS credentials](https://docs.aws.amazon.com/cli/latest/userguide/cli-configure-files.html) files.

Examples of both are provided in the variant_caller repository at `./serverless_genomics/variant_caller/varcall_client`.

To upload them, you can simply use `scp`, pointing to your private key file (the same used for SSH access to the VM):
```
#upload aws credentials
scp -i private_key.pem <1: 1 windows (created Wed Jun 29 17:31:54 2022)
/path/to/.aws/credentials> ubuntu@<VirtualMachine_DNS>:/root/aws_credentials

# upload lithops config
scp -i private_key.pem </path/to/.lithops/config> ubuntu@<VirtualMachine_DNS>:/root/lithops_config
```

### 2.4 Running the Docker container
From the **VM** terminal, you can then run the Docker. We have provided a convenience script for you, which "mounts" the necessary files from the EC2 environment so that they are accessible within the container.

These are:

| EC2 Path                         | Container Path            | Description                                                                                                                                                          |
|----------------------------------|---------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| /tmp                             | /tmp                      | This is where the lithops logs will be written (to /tmp/lithops/logs)                                                                                                |
| /home/ubuntu/serverless_genomics | /root/serverless_genomics | Copy of the variant_caller gitlab repository. This contains the code for running the Docker container, as well as the code for running the pipeline                  |
| /home/ubuntu/aws_credentials     | /root/credentials         | AWS credentials file (containing your super-secret key information). Inside the container, this is set to the environmental variable AWS_SHARED_CREDENTIALS_FILE     |
| /home/ubuntu/lithops_config      | /root/config              | Lithops configuration file, containing the details of your serverless execution. Inside the container, this is set to the environmental variable LITHOPS_CONFIG_FILE |

You can run the Docker from the **VM** home directory simply by providing these paths as arguments for the convenience script (note you also provide the container image name, e.g., tkchafin/varcall_client:0.2):
```
./serverless_genomics/variant_caller/varcall_client/varcall_docker_run.sh ./serverless_genomics/ ./aws_config ./litops_config tkchafin/varcall_client:0.2
```

That's it! You now should notice your terminal prompt has changed to something like:
```
root@4d311b7ae8e8:~/serverless_genomics#
```

You can now run the pipeline.

###2.4b Running the Docker container in a background shell

Although you can now run the pipeline, if you logout or exit your **VM** terminal, any ongoing pipeline runs or running images will stop. To prevent this, you can run a background shell using the tool `tmux`:

```
# start a tmux instance
# should notice a green bar at bottom of terminal after running
tmux

# run the docker
./serverless_genomics/variant_caller/varcall_client/varcall_docker_run.sh ./serverless_genomics/ ./aws_config ./litops_config tkchafin/varcall_client:0.2
```

As above, your terminal prompt should change to indicate you are "in" the container. To detach from the `tmux` instance, just type `CTRL+B`, then `D`. Now, you can safely step away for a coffee and return to have mapped and variant-called data!

To revisit your running tmux sessions, you can check for running instances using `tmux ls`:
```
ubuntu@ip-172-31-82-75:~$ tmux ls
0: 1 windows (created Wed Jun 29 17:18:28 2022)
ubuntu@ip-172-31-82-75:~$
```

Now, "re-attach" using `tmux attach`:
```
tmux attach -t 0
```

## 3 Building the container from scratch

If you have an issue with Docker (and can't wait for help by posting over at the [Gitlab Issues page](https://gitlab.bioss.ac.uk/lmarcello/serverless_genomics/-/issues)), or want to add some functionality to the Docker, you can find the complete Dockerfile at `serverless_genomics/variant_caller/varcall_client/Dockerfile`. Inside, you will find directions for building the container.

To build the docker from scratch (and with docker installed on your **local** machine), you can type:
```
# e.g. docker build -t tkchafin/varcall_client:0.2 varcall_client/
docker build -t <username>/varcall_client:<tag> varcall_client/
```
Note that this can take a while, especially if it is the first time you are building it (hence no cache for docker to use). 
