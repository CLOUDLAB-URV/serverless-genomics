#!/bin/bash
# wrapper script to converts a .map alignment file to mpileup format

map_file="${1:-/dev/stdin}"
ref_genome="${2:-/dev/stdin}"
#echo "Starting gempileup_run.sh"
bash /function/bin/gempileup_v7.sh $map_file $ref_genome   | 
bash  /function/bin/gempileup_merge.sh > $map_file".mpileup"
#echo "End of gempileup_run.sh"