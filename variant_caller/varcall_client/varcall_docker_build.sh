#!/bin/bash

WORKDIR=/home/tchafin/bioss/serverless_genomics/variant_caller

cd $WORKDIR

docker build -t tkchafin/varcall_client:0.2 varcall_client/
