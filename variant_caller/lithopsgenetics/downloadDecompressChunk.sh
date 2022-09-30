#!/bin/bash
# Script downloadDecompressChunks.sh: Dowloads a byte range chunk from a gzip file stored in a BUCKET
#                                     unzips it and stores it back in the BUCKET

# COMMAND LINE INPUT
file=$1
link=$2
bucket=$3
start_line=$4
block_length=$5
start_byte=$6
end_byte=$7

# 1. DOWLOAD THE CHUNK FROM THE BUCKET
start_byte_=$(($start_byte-1))
curl -r "${start_byte_}-${end_byte}" "${link}${file}" > "${file}_${start_line}.fastq.gz"

# 2. UNZIP THE CHUNK
gztool -I "${file}i" -n ${start_byte} -L ${start_line} "${file}_${start_line}.fastq.gz" | head -${block_length} > "${file}_${start_line}.fastq"


# 3. UPLOAD FILES TO THE BUCKET
lithops storage put "${file}_${start_line}.fastq" $bucket

# 4. DELETE CHUNK TO FREE UP SPACE
rm "${file}_${start_line}.fastq.gz"
rm "${file}_${start_line}.fastq"
