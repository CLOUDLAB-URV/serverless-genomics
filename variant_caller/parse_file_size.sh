#!/bin/bash
# parse multiple logs to compare specific file sizes across runs
cd varcall_out || exit
# inputs
subset=""  # subset size for fq or fa input file, parsed from log name
input_type="fq"  # fq or fa
target_file="_filt_wline_no_corrected.map.mpileup.csv"

target_file_type="${target_file##*.}"
if [[ $input_type = "fq" ]]; then
    grep_target=${subset}"Kfq"
fi
if [[ $input_type = "fa" ]]; then
    grep_target=${subset}"Mfa"
fi
 
log_subset=*${grep_target}*rep1.log

#example
#grep _filt_wline_no_corrected.map.mpileup.csv  *1000Kfq*rep1.log |
grep $target_file $log_subset |
awk -v target_file_type="$target_file_type" '
{
    split($1, log_name_parts, /_/)

    # fastq size
    fastq_size=log_name_parts[5]
    fastq_size=gensub(/Kfq/, "","g",fastq_size)
    
    # fasta size
    fasta_size=log_name_parts[6]
    fasta_size=gensub(/Mfa/, "","g",fasta_size)
    
    print $1"\t"fastq_size"\t"fasta_size"\t"target_file_type"\t"$7

}
' > ${target_file_type}"_sizes_"${subset}${input_type}"_subset.txt"