#!/bin/bash
# program that filters unmatched indices (reads with alignment(s) only in one chunk), to return only reads with matches in multiple fasta chunks
file=$1
split=$2
echo "processing "$file
out_name=$split

awk -v OUT="$out_name" -v DIR="$1" '
{
	## SELECT INDICES WITH A MATCH

	if ($3=="2F") {
                printf ("%s\t%s\n", $1, $2) > OUT
        }
}
' $1

if [[ ! -f $out_name ]]
then
    touch $out_name
fi 