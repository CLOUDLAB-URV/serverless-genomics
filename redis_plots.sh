#!/bin/bash
# script to make plots for the redis log
log=$1
dir="varcall_out/"
dir=""
start_log=${dir}${log}"_redis_start_times.txt"
failed_functions=${dir}${log}"_redis_failed_functions.txt"
exec_times=${dir}${log}"_redis_execution_time.txt"

# start time table
grep "redis: time_start"  ${dir}${log} | awk '{time=gensub("time_start:", "","g", $4); print $2"\t"time}' > $start_log

# get functions that fail
grep "corrected map index filename:" ${dir}${log} | awk '{print $2}' | sort -k1n | awk '{for(i=p+1; i<$1; i++) print i} {p=$1}'  >  $failed_functions

# execution tim
grep "redis: execution_time:" ${dir}${log} | awk '{time=gensub(":", "","g", $5); print $2"\t"time}' > $exec_times

# Run R script
TZ="Australia/Sydney" Rscript "redis_plots.R" $start_log $failed_functions $exec_times
