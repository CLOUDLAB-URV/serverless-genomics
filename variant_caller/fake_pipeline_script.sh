#!/bin/bash

file="varcall_out/varcall_1_1_SRR15068323_800Kfq_100Mfa_5itn_rep1.log"
while read line; do
# reading each line
echo $line
done < $file