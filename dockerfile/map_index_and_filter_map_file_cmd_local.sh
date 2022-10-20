#!/bin/bash
# script that takes stdout from gem mapper (with map flag) and
#/bin/gem3-mapper/bin/gem-mapper -I /home/lumar/bioss/test_data/Hsapiens/fasta/hg19_145450269split_00000001.gem -i /home/lumar/bioss/test_data/Hsapiens/fastq/split_22/SRR15068323_split_1.fastq  -F MAP | bash /home/lumar/bioss/cloudbutton/scripts/07_aln_score_parser/parse_gem_maxindex_minimapfile_stdin.sh  /dev/stdin SRR15068323_split_1_22__hg19_145450269split_00000001.se
#map_index_and_filter_map_file_cmd.sh /home/lumar/bioss/test_data/Tbrucei/fasta/TriTrypDB-9.0_TbruceiTREU927_Genome_MbChr_split_1.gem /home/lumar/bioss/test_data/Tbrucei/fastq/SRR6052133_1_split_1.fastq /home/lumar/bioss/test_data/Tbrucei/fastq/SRR6052133_2_split_1.fastq SRR6052133
gem_ref=$1
fastq=$2
fastq2=$3
map_file_nosuffix=$4
seq_type=$5
echo "gem-mapper / aln indexer program starting"
echo "gem reference: "$gem_ref
echo "fastq file: "$fastq
echo "fastq2 file: "$fastq2
echo "map file no suffix: "$map_file_nosuffix


if [ $seq_type == "single-end" ]
then
    echo "single-end alignment" 
    gem-mapper -I "$gem_ref" -i "$fastq" -F MAP | bash /function/bin/parse_gem_maxindex_minimapfile_stdin_v2.sh  /dev/stdin $map_file_nosuffix $seq_type
else
    echo "paired-end alignment"
    gem-mapper -I "$gem_ref" -1 "$fastq" -2 "$fastq2"  -F MAP | bash /function/bin/parse_gem_maxindex_minimapfile_stdin_v2.sh  /dev/stdin $map_file_nosuffix $seq_type
fi



