#!/bin/bash

if [ "$#" -ne 4 ]; then
  echo "Usage: $0 </path/to/serverless_genomics> </path/to/.aws/credentials> </path/to/.lithops/config> <DOCKER_IMAGE>" >&2
  exit 1
fi
if ! [ -e "$2" ]; then
  echo "$2 not found" >&2
  echo "Usage: $0 </path/to/serverless_genomics> </path/to/.aws/credentials> </path/to/.lithops/config>" >&2
  exit 1
fi
if ! [ -d "$1" ]; then
  echo "Usage: $0 </path/to/serverless_genomics> </path/to/.aws/credentials> </path/to/.lithops/config>" >&2
  echo "$1 not a directory" >&2
  exit 1
fi
if ! [ -e "$3" ]; then
  echo "$3 not found" >&2
  echo "Usage: $0 </path/to/serverless_genomics> </path/to/.aws/credentials> </path/to/.lithops/config>" >&2
  exit 1
fi

#set manually
#Change me!
#WORKDIR=/home/tchafin/bioss/serverless_genomics
#MOUNT=/root/bioss/serverless_genomics
#AWS=/home/tchafin/.aws/credentials
#LITHOPS=/home/tchafin/.lithops/config

WORKDIR=$(realpath "$1")
AWS=$(realpath "$2")
LITHOPS=$(realpath "$3")
DOCKER=$4

cd $WORKDIR
MOUNT=/root/serverless_genomics
docker run -it -v $WORKDIR:$MOUNT \
    -v /tmp:/tmp \
    -w $MOUNT \
    -v $AWS:/root/.aws/credentials \
    -v $LITHOPS:/root/.lithops/config \
    $DOCKER
