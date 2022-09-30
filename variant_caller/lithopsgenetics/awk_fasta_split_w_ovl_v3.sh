#!/bin/bash
# script that takes a fasta file and splits it into smaller files

file=$1
block_length=$2
ovl=$3
file_name="${file%.*}"

awk -v BLOCK_LEN="$block_length" -v OVERLAP="$ovl" -v NAME="$file_name" ' 
function process_sequence()
    {
        while ((l=length(seq))>0) 
            {
                if (l>=left) {
                    print ">"name"_"(++cntr)"\n"substr(seq,1,left+OVERLAP) > file; 
                    close(file); 
                    file=NAME"_split_"sprintf("%08d",++chunks)".fasta"; 
                    seq=substr(seq,left+1); 
                    left=BLOCK_LEN
                } 
                else {
                    print ">"name"_"(++cntr)"\n"seq > file; 
                    seq=""; 
                    left=left-l
                }
            }
        } 
BEGIN{
    file=NAME"_split_"sprintf("%08d",++chunks)".fasta"; 
    seq=""; 
    left=BLOCK_LEN} 
    {
        if ($0~"^>") {
            process_sequence(); 
            name=substr($0,2)
        } 
        else {
            seq=seq $0
        }
    } 
END{process_sequence()}' $file 
