#!/bin/bash
# script that wraps gem mapper and parse_gem_maxindex_minimapfile_stdin_v2.sh to allow the command to work using subprocess.
#/bin/gem3-mapper/bin/gem-mapper -I /home/lumar/bioss/test_data/Hsapiens/fasta/hg19_145450269split_00000001.gem -i /home/lumar/bioss/test_data/Hsapiens/fastq/split_22/SRR15068323_split_1.fastq  -F MAP | bash /home/lumar/bioss/cloudbutton/scripts/07_aln_score_parser/parse_gem_maxindex_minimapfile_stdin_v2.sh  /dev/stdin SRR15068323_split_1_22__hg19_145450269split_00000001.se
gem_ref=$1
fastq=$2
fastq2=$3
map_file_nosuffix=$4
database=$5
seq_type=$6
threads=$7
echo "gem-mapper / aln indexer program starting"
echo "gem reference: "$gem_ref
echo "fastq file: "$fastq
echo "fastq2 file: "$fastq2
echo "map file no suffix: "$map_file_nosuffix
if [ "$database" == "s3" ]
then
    echo "database: "$database
    if [ $seq_type == "single-end" ]
    then
        echo "single-end alignment - file from "$database
        gem-mapper -I "$gem_ref" -i "$fastq" -F MAP -t "$threads" | bash /function/bin/parse_gem_maxindex_minimapfile_stdin_v2.sh  /dev/stdin $map_file_nosuffix $seq_type
    else
        echo "paired-end alignment - file from "$database
        gem-mapper -I "$gem_ref" -1 "$fastq" -2 "$fastq2" -F MAP -t "$threads" | bash /function/bin/parse_gem_maxindex_minimapfile_stdin_v2.sh  /dev/stdin $map_file_nosuffix $seq_type
    fi
elif  [ "$database" == "SRA" ]
then 
    echo "database: "$database
    if [ $seq_type == "single-end" ]
    then
        echo "single-end alignment - file from "$database
        gem-mapper -I "$gem_ref" -i "$fastq" -F MAP -t "$threads" | bash /function/bin/parse_gem_maxindex_minimapfile_stdin_v2.sh  /dev/stdin $map_file_nosuffix $seq_type
    else
        echo "paired-end alignment - file from "$database
        gem-mapper -I "$gem_ref" -i "$fastq" -p -F MAP -t "$threads" | bash /function/bin/parse_gem_maxindex_minimapfile_stdin_v2.sh  /dev/stdin $map_file_nosuffix $seq_type

    fi
fi
