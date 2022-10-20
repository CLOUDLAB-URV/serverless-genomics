#!/bin/bash
# Script generateUploadIndexInfoBucket.sh: Takes a gzip file and creates an index for it, reformates the index info,
#                                          stores important files in a bucket

# COMMAND LINE INPUT
file=$1

# 1. INDEX FASTQ.GZ AND REFORMAT INDEX INFO
# Create index -i, with line number information -x, with span in uncompressed MiB between index points
gztool -i -x -s 1 $file
# Produce index info: showing internals of all index files in this directory -e, showing data about each index point
gztool -ell "${file}i" > "${file}i.info"
# Get total lines and size form filei.info
tot_lines=$(awk ' /Number of lines/ { print $5 } ' "${file}i.info")
tot_size=$(awk ' /Guessed gzip/ { print $10 } ' "${file}i.info" | awk '{gsub(/^.{1}/,"");}1')
# add 1 to tot_size, as in all blocks 1 is subtracted to match end of previous block
tot_size="$((tot_size + 1))"
# reformat .info file to have all points in list start at new line and remove header
sed -e 's/#/\n#/g' -e 's/L//g'  "${file}i.info" |
    awk ' /^#[0-9]+/ { printf ("%s %s\n", $3, $6);} ' > "${file}i_tab.info"
echo $tot_size" "$tot_lines >> "${file}i_tab.info"

echo "Total lines:"$tot_lines
echo "Total size:"$tot_size
