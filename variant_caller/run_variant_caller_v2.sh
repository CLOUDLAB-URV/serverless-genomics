#!/bin/bash
# script to run the variant caller

function join_by {
  local d=${1-} f=${2-}
  if shift 2; then
    printf %s "$f" "${@/#/$d}"
  fi
}

input=$1  #file containing command line cmd_args and associated values as tab separated file with args on one line and values on another (in columns)
run=$2
loop=$3
run_script=$4


echo "run: "$run
echo "loop: "$loop
# running script or only wrapper
echo "running the script: "$run_script

out_dir="varcall_out/"
#create output directory if not present
mkdir -p $out_dir

# make sure sratoolkit path is accessible
export PATH="$PATH:/bin/sratoolkit.3.0.0-ubuntu64/bin"
# remove carriage return from input (if open in excel)
awk '{print gensub("\r$","",1)}' $input > ${input}"_awk"
mv ${input}"_awk" $input
# create command_line
cmd_args=$(cat $input | grep -v ^\# | head -1 | tail -1)
cmd_vals=$(cat $input | grep -v ^\# | head -2 | tail -1)
# convert to array
IFS=$'\t' read -r -a arg_array <<< "$cmd_args"
IFS=$'\t' read -r -a val_array <<< "$cmd_vals"
test_n=${val_array[0]}
command="python varcall_lithops_demo_v7.py"
for ((i=1;i<${#arg_array[@]};i++)); 
do 
    # generate command line
    if [ "${val_array[$i]}" = "NA" ]; then
        #echo "argument not in use:  -"${arg_array[$i]}" "${val_array[$i]}
        continue
    else
        command=$command" -"${arg_array[$i]}" "${val_array[$i]}
        echo " -"${arg_array[$i]}" "${val_array[$i]}
    fi

    # generate components for output name
    if [ "${arg_array[$i]}" = "fq" ]; then
        fq=${val_array[$i]}
    elif [ "${arg_array[$i]}" = "nfq" ]; then
        nfq=${val_array[$i]}
        nfq=$(echo $(( nfq / 1000 )))"Kfq"
    elif [ "${arg_array[$i]}" = "nfa" ]; then
        nfa=${val_array[$i]}
        nfa=$(echo $(( nfa / 1000000 )))"Mfa"
    elif [ "${arg_array[$i]}" = "itn" ]; then
        itn=${val_array[$i]}"itn"
    fi
done

echo "command "$command
echo "fq: "$fq
echo "nfq: "$nfq
echo "nfa: "$nfa
echo "itn: "$itn
out_name="varcall_"$test_n"_"$run"_"$fq"_"$nfq"_"$nfa"_"$itn
echo "output file name: "$out_name


# create log files (across all iterations)
exec_time_log=$out_dir${out_name}"_exec_time.log"
file_size_log=$out_dir${out_name}"_file_size.log"
>$exec_time_log
>$file_size_log
 
echo -e "running $loop command iterations" 
for ((i=1;i<=$loop;i++)); 
do 
    ############### 1. RUN SCRIPT
    echo -e "\n\n STARTING ITERATION "$i
    # execute script and write general log to output directory (one file per iteration)
    # if debugging a log file don't run the script, just read the log
    # run script, combine stdout and stdout print to file and to screen
    
    if [ $run_script = "True" ]
    then
        echo "running the variant caller" 
        iter_log=$out_dir${out_name}"_rep"${i}".log"
        graphics_log=$out_dir${out_name}"_rep"${i}"_graphics.log"
        unbuffer $command 2>&1 | tee $iter_log | ./graphics_log_input_creator > $graphics_log
        
    else
        echo "variant caller log parsing"
        iter_log=$out_dir${out_name}"_rep"${i}".log"
    fi
    
    echo -e "\n\n END OF ITERATION "$i
    
    if test -f "$iter_log"; then
        # ############### 2. REFORMAT LOG
        # echo "reformatting log"

        # iter_log_sorted="${iter_log}_functions_sorted.log"
        # >$iter_log_sorted
        
        # awk '
        # { # adding 4 digit random number
        #     if ($0~"^[A-Z][0-9]{4}\t") {
        #             what=gensub("^([A-Z][0-9]{4}\t[0-9]+).*$","\\1",1) 
        #             message=gensub("^[A-Z][0-9]{4}\t[0-9]+\t","",1) 
        #             if (message!="") t[what]=t[what]message"\n"
        #     } else {
        #         if ($0!="") print
        #     }
        # } 
        # # {  # previous version - before adding 4 digit random number
        # #     if ($0~"^[A-Z]\t") {
        # #             what=gensub("^([A-Z]\t[0-9]+).*$","\\1",1) 
        # #             message=gensub("^[A-Z]\t[0-9]+\t","",1) 
        # #             if (message!="") t[what]=t[what]message"\n"
        # #     } else {
        # #         if ($0!="") print
        # #     }
        # # } 
        # END {
        #     for (i in t) {
        #         l=split(t[i],s,"\n") 
        #         for (j=1;j<=l;++j) {
        #             print i"\t"s[j]
        #         }
        #     }
        # }
        # ' $iter_log > $iter_log_sorted

        iter_log_sorted=$iter_log
        
        ############### 3. EXTRACT TIME AND SIZE DATA
        echo "extracting execution times from log files"
        grep "execution_time" $iter_log_sorted | sed "s/^/$i\t/" | awk '{out=gensub(":","\t","g",$0); print out}'>>  $exec_time_log
        grep "time_start" $iter_log_sorted | sed "s/^/$i\t/" | awk '{out=gensub(":","\t","g",$0); print out}'>>  $exec_time_log
        grep "time_end" $iter_log_sorted | sed "s/^/$i\t/" | awk '{out=gensub(":","\t","g",$0); print out}'>>  $exec_time_log


        echo "extracting file sizes from log files"
        awk -v iter_n="$i" 'NF==11 && match($11, /fastq$/) {if (!a[$0]++) {print iter_n"\t"$1"\t"$2"\t"substr($11, RSTART, RLENGTH)"\t"$11"\t"$7}}' $iter_log  >>  $file_size_log
        awk -v iter_n="$i" 'NF==11 && match($11, /fasta$/) {if (!a[$0]++) {print iter_n"\t"$1"\t"$2"\t"substr($11, RSTART, RLENGTH)"\t"$11"\t"$7}}' $iter_log  >>  $file_size_log
        awk -v iter_n="$i" 'NF==11 && match($11, /map.index.txt$$/) {if (!a[$0]++) {print iter_n"\t"$1"\t"$2"\t"substr($11, RSTART, RLENGTH-4)"\t"$11"\t"$7}}' $iter_log  >>  $file_size_log
        awk -v iter_n="$i" 'NF==11 && match($11, /map.corrected_index.txt$/) {if (!a[$0]++) {print iter_n"\t"$1"\t"$2"\t"substr($11, RSTART, RLENGTH-4)"\t"$11"\t"$7}}' $iter_log  >>  $file_size_log
        awk -v iter_n="$i" 'NF==11 && match($11, /filt_wline_no.map$/) {if (!a[$0]++) {print iter_n"\t"$1"\t"$2"\t"substr($11, RSTART+14, RLENGTH)"\t"$11"\t"$7}}' $iter_log  >>  $file_size_log
        awk -v iter_n="$i" 'NF==11 && match($11, /filt_wline_no_corrected.map$/) {if (!a[$0]++) {print iter_n"\t"$1"\t"$2"\t"substr($11, RSTART+14, RLENGTH)"\t"$11"\t"$7}}' $iter_log  >>  $file_size_log
        awk -v iter_n="$i" 'NF==11 && match($11, /filt_wline_no_corrected.map.mpileup$/) {if (!a[$0]++) {print iter_n"\t"$1"\t"$2"\t"substr($11, RSTART+28, RLENGTH)"\t"$11"\t"$7}}' $iter_log  >>  $file_size_log
        awk -v iter_n="$i" 'NF==11 && match($11, /filt_wline_no_corrected.map.mpileup.csv$/) {if (!a[$0]++) {print iter_n"\t"$1"\t"$2"\t"substr($11, RSTART+28, RLENGTH)"\t"$11"\t"$7}}' $iter_log   >>  $file_size_log
        awk -v iter_n="$i" 'NF==11 && match($11, /filt_wline_no_corrected.map.mpileup.parquet$/) {if (!a[$0]++) {print iter_n"\t"$1"\t"$2"\t"substr($11, RSTART+28, RLENGTH)"\t"$11"\t"$7}}' $iter_log   >>  $file_size_log
        awk -v iter_n="$i" 'NF==11 && match($11, /reduce.mpileup$/) {if (!a[$0]++) {print iter_n"\t"$1"\t"$2"\t"substr($11, RSTART, RLENGTH)"\t"$11"\t"$7}}' $iter_log   >>  $file_size_log
        awk -v iter_n="$i" 'NF==11 && match($11, /merged.mpileup$/) {if (!a[$0]++) {print iter_n"\t"$1"\t"$2"\t"substr($11, RSTART, RLENGTH)"\t"$11"\t"$7}}' $iter_log   >>  $file_size_log

        ############### 4. GENERATE SUMMARY STATS
        echo -e "\n#######END OF LOG\n"
        echo -e "\ncommand line used: \n"$command >> $iter_log_sorted
        echo -e "\nSUMMARY STATS" >> $iter_log_sorted
        echo -e "\nMAP STAGE" >> $iter_log_sorted
        iterdata_n=$(grep "number of iterdata elements" $iter_log_sorted | awk '{print $5}') 
        echo -e "Functions launched:\t"$iterdata_n >> $iter_log_sorted
        echo -e "Functions per step" >> $iter_log_sorted

        gem_ind_yes=$(grep "was successfully built" $iter_log_sorted | awk '!arr[$2]++ {print $2}' | wc -l)
        echo -e "1. gem indexer:\t"$gem_ind_yes >> $iter_log_sorted
        
        mapper_yes=$(awk 'NF==11 && match($11, /filt_wline_no.map$/) {if (!a[$0]++) {print $0}}' $iter_log_sorted | sort --uniq | wc -l)
        echo -e "2. gem mapper:\t"$mapper_yes >> $iter_log_sorted

        redis_found=$(grep "found corrected index" $iter_log_sorted | sort --uniq | wc -l) 
        echo -e "3a. redis index found:\t"$redis_found >> $iter_log_sorted

        redis_download=$(grep "corrected map index filename:" $iter_log_sorted | sort --uniq | wc -l) 
        echo -e "3b. redis index downloaded:\t"$redis_download >> $iter_log_sorted

        map_corrected=$(grep "_filt_wline_no_corrected.map$" $iter_log_sorted | sort --uniq | wc -l) 
        echo -e "4. correct map file:\t"$redis_download >> $iter_log_sorted

        mpileup_yes=$(grep "_filt_wline_no_corrected.map.mpileup$" $iter_log_sorted | sort --uniq | wc -l) 
        echo -e "5. gempileup:\t"$mpileup_yes >> $iter_log_sorted

        csv_yes=$(grep "filt_wline_no_corrected.map.mpileup.csv" $iter_log_sorted | sort --uniq| wc -l)
        echo -e "6. csv:\t"$csv_yes >> $iter_log_sorted
        
        echo -e "\nREDUCE (STAGE 1)" >> $iter_log_sorted
        reduce1_n=$(grep "reduce function 1" $iter_log_sorted | sort --uniq | wc -l)
        echo -e "Functions launched:\t"$reduce1_n >> $iter_log_sorted
        
        echo -e "\nREDUCE (STAGE 2)" >> $iter_log_sorted
        reduce2_n=$(grep "reduce function 2" $iter_log_sorted | sort --uniq | wc -l)
        echo -e "Functions launched:\t"$reduce1_n >> $iter_log_sorted
        echo -e "Functions per step" >> $iter_log_sorted
 
        reduce2_input=$(grep "reduce.mpileup$" $iter_log_sorted | sort --uniq | wc -l)
        echo -e "1. reduced mpileup input:\t"$reduce2_input >> $iter_log_sorted

        reduce2_output=$(grep "reduce.mpileup_merged.mpileup$" $iter_log_sorted | sort --uniq | wc -l)
        echo -e "2. merged mpileup output:\t"$reduce2_output >> $iter_log_sorted

        #####
        iter_log_summary="varcall_results_summary.tsv" # for all varcall runs
        #iter_log_summary_header=$(head -n 1 $iter_log_summary)
        repeat_n=$i
        out_args=$(join_by '\t' $cmd_args "repeat_n" "main_log"      "iterdata_n" "gem_ind_yes" "mapper_yes" "redis_found" "redis_download" "map_corrected" "mpileup_yes" "csv_yes" "reduce1_n" "reduce2_n" "reduce2_input" "reduce2_output" "command")
        out_vals=$(join_by '\t' $cmd_vals $repeat_n $iter_log_sorted $iterdata_n  $gem_ind_yes  $mapper_yes  $redis_found   $redis_download $map_corrected  $mpileup_yes  $csv_yes  $reduce1_n  $reduce2_n  $reduce2_input  $reduce2_output  "$command")
        if  [[ -z $(grep '[^[:space:]]' $iter_log_summary) ]] ; then # if file is empty print header
            echo -e $out_args >> $iter_log_summary
        fi
        echo -e $out_vals >> $iter_log_summary
        grep DEBUG_TEST  $iter_log_sorted
        grep execution_time_total_varcall  $iter_log_sorted
        tail -n 22 $iter_log_sorted

        bash redis_plots.sh  $iter_log_sorted
    fi
done

echo -e "\niterations complete" grep execution_time_total_varcall  $iter_log_sorted


echo -e "\n Running R script process_execution_times.R for execution time parsing"
command_for_ggplot="$(awk -v cmd="$command" 'BEGIN{cmd = gensub(/-fa/,"\n-fa","g",cmd); cmd = gensub(/-nlfq/,"\n-nlfq","g",cmd); print cmd }')"
echo -e "command for ggplot \n""$command_for_ggplot"

TZ="Australia/Sydney" Rscript $prog_dir"process_execution_times.R" $exec_time_log $file_size_log "$command_for_ggplot"
