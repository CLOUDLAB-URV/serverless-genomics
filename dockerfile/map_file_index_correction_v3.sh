#!/bin/bash
# script that takes as input:
# 1. filtered map file (with original line numbers, and best index field) as obtained from parse_gem_maxindex_minimapfile_stdin.sh
# 2. The corrected index file
# and removes any alignment that is below the threshold (max index value + tolerance) in terms of its index value.
# threshold is defined as "max index + tolerance" in the case of reads with no matches elsewhere, 
# or as "max corrected index + tolerance" for reads with matches in different chunks

# command line example 
#bash map_file_index_correction_v3.sh  /home/lumar/bioss/cloudbutton/scripts/SRR6052133_split_12_filt_wline_no.map /home/lumar/bioss/cloudbutton/scripts/SRR6052133_split_corrected_index.txt paired-end 0 False

#echo "map_file_index_correction_v3.sh starting"

# ARG 1: filtered map file 

map_file_default="/home/lumar/bioss/test_data/Hsapiens/map/SRR15068323_split_1_22__hg19_145450269split_00000001.se_filt_wline_no.map"
map_file=${1:-$map_file_default}

map_file_name="${map_file%.*}"

# ARG 2: corrected index file
index_file_default="/home/lumar/bioss/cloudbutton/scripts/07_aln_score_parser/merged_map_1_2_3_4_shared_only.txt"
index_file=${2:-$index_file_default}

# ARG 3: paired- v single-end sequencing
seq_type=${3:-"single-end"}

# ARG 4: tolerance value
# specify how many strata to take below top index
tolerance=${4:-0}

# ARG 5: debug option
debug="${5:-"False"}"

echo -e "processing $seq_type map file: $map_file\nusing corrected index file: $index_file\nwith a tolerance of **$tolerance** additional strata."

awk -v NAME="$map_file_name" -v seq_type="$seq_type" -v TOLERANCE="$tolerance" -v debug="$debug" '
function print_debug(debug, string)
{
    if (debug =="True") {
        print string
    }
}
function printf_debug(debug, format, var1, var2, var3)
{
    if (debug =="True") {
        printf (format, var1, var2, var3)
    }
}
function join(array, start, end, sep,    result, i)  #https://www.gnu.org/software/gawk/manual/html_node/Join-Function.html
{
    if (sep == "")
        sep = " "
    else if (sep == SUBSEP) # magic value
        sep = ""
    result = array[start]
    for (i = start + 1; i <= end; i++)
        result = result sep array[i]
    return result
}
function count_ok_aln(match_summary_pt1, threshold)
{
    stratum_no = split(match_summary_pt1, strata, /:/)
    # number of alignments to keep (above or equal to threshold)
    ok_aln = 0
    for (i = 1; i <= stratum_no; ++i) {
        print_debug(debug,"index position: "i" index value: "strata[i])
        if (strata[i]~"[1-9]" && i <= threshold) {
                print_debug(debug, "found alignment above threshold")
                ok_aln = ok_aln + strata[i]     #0:0:1:1:0  - define maximum index as $2 + tolerance
                print_debug(debug, "interim number of ok alignments: "ok_aln)
        }
    }
    return ok_aln
}
function filter_aln(alignments, ok_aln)
{
    align_data = split(alignments, aligns, /,/)
    for (i = 1; i <= align_data; ++i) {
        if (i > ok_aln) {
            delete aligns[i]
        }
    }
    filtered_aligns = join(aligns, 1, length(aligns),",")
    return filtered_aligns
}
BEGIN {
    #counts
    count_reads=0 # all reads
    count_match=0 # all matching reads
    count_filter=0 # matching reads that are filtered
    count_unchanged=0 # matching reads that are unchanged
    count_discard=0  # matching reads that are discarded
    count_other=0 # non matching reads
    count_other_filter=0 # other reads that are filtered
    count_other_unchanged=0 # other reads that are unchanged
}
{ 
if (NR==FNR) {
        t[$1]=$0
}
else {
    split(t[$1],s,"\t");
    count_reads++
    if (count_reads == 1) {
        print_debug(debug,"column 1: "$1)
        print_debug(debug,"column 2: "$2)
        print_debug(debug,"column 3: "$3)
        print_debug(debug,"column 4: "$4)
        print_debug(debug,"column 5: "$5)
        print_debug(debug,"column 6: "$6)
        print_debug(debug,"column 7: "$7)
        print_debug(debug,"column 8: "$8)
        print_debug(debug,"column 9: "$9)
    }
    
    #print_debug(debug, "number of fields: "NF)
    if (seq_type == "paired-end" && NF<9) {
        #print_debug(debug, "read is not paired end - discarded")
        next
    }
    read_name=$3
    if (seq_type == "single-end") {
        read=$4
        qual=$5
        summary=$6
        aln_col=$7
    }
    else { # paired-end
        read1=$4
        read2=$5
        qual1=$6
        qual2=$7
        read=read1" "read2
        qual=qual1" "qual2
        summary=$8
        aln_col=$9
    }
    # match summary field
      
    # split match summary into two parts (before and after +)
    split(summary, summary_parts, /\+/)
    match_status = ""

    #print_debug(debug,"index matching?: "$1"\t"s[1])
    # if line number matches, compare max corrected index (corrected index - s[2]) with max read index ($2) in chunk
    if ($1 == s[1]) { 
        print_debug(debug, "\nread no. "count_reads)
        print_debug(debug, $0)
        print_debug(debug,"match summary: "summary)
        print_debug(debug,"match summary part 1: "summary_parts[1])
        print_debug(debug, "alignments: "aln_col)
        print_debug(debug,"found match: "$1"\t"s[1])
        count_match++
        threshold = s[2] + TOLERANCE
        print_debug(debug, "threshold: "threshold", alignment index: "$2)
        if ($2 <= threshold) {
            match_status="keep"
            # find number of alignments to keep based on threshold
            ok_aln = count_ok_aln(summary_parts[1], threshold)
            print_debug(debug,"number of ok alignments (from summary): "ok_aln)
            # find number of alignments listed in alignment field
            tot_aln = split(aln_col, aligns, /,/)
            print_debug(debug,"number of total alignments (from aln coln): "tot_aln)
            # if total alignments to keep < total number of alignments, then remove the unwanted alignments.
            if (tot_aln > ok_aln) {
                print_debug(debug, "match - filtered alignments")
                count_filter++
                filtered_aligns=filter_aln(aln_col, ok_aln)
                print_debug(debug, t[$1]"\t"$1"\t"$2"\t"match_status"\t"summary"\t"ok_aln"\t"tot_aln"\t"aln_col)
                #print t[$1]"\t"$1"\t"$2"\t"match_status"\t"summary"\t"ok_aln"\t"tot_aln"\t"filtered_aligns
                printf ("%s\t%s\t%s\t%s\t%s\n",read_name, read, qual, summary, filtered_aligns)  > NAME"_corrected.map"
                print_debug(debug, read_name"\t"read"\t"qual"\t"summary"\t"filtered_aligns)
            }
            else {
                print_debug(debug, "match - alignments unchanged")
                count_unchanged++
                printf ("%s\t%s\t%s\t%s\t%s\n", read_name, read, qual, summary, aln_col)  > NAME"_corrected.map"
                print_debug(debug, read_name"\t"read"\t"qual"\t"summary"\t"aln_col)
            }
        }
        else {
            # discard read from chunk
            count_discard++
            match_status = "discard"
            print_debug(debug,"match - read discarded - no alns at given threshold: "t[$1]"\t"$1"\t"$2"\t"match_status"\t"summary)
            print_debug(debug, read_name"\t"read"\t"qual"\t"summary"\t"aln_col)

        }
    }
    else {
        
        # if line number does not match, output as .map file omitting line number and index, and use threshold = "index + tolerance" 
        count_other++
        threshold = $2 + TOLERANCE
        ok_aln = count_ok_aln(summary_parts[1], threshold)
        # find number of alignments listed in alignment field
        tot_aln = split(aln_col, aligns, /,/)
        # if total alignments to keep < total number of alignments, then remove the unwanted alignments.
        if (tot_aln > ok_aln) {
            print_debug(debug, "no match - filtered alignments")
            count_other_filter++
            filtered_aligns=filter_aln(aln_col, ok_aln)
            print_debug(debug, "filtered alignments: "filtered_aligns)
            printf ("%s\t%s\t%s\t%s\t%s\n", read_name, read, qual, summary, filtered_aligns)  > NAME"_corrected.map"
            print_debug(debug, read_name"\t"read"\t"qual"\t"summary"\t"filtered_aligns)
        }
        else {
            print_debug(debug, "no match - unchanged alignments")
            count_other_unchanged++
            printf ("%s\t%s\t%s\t%s\t%s\n",read_name, read, qual, summary, aln_col)  > NAME"_corrected.map"
            print_debug(debug, read_name"\t"read"\t"qual"\t"summary"\t"aln_col)
        }

    } 
}
}
END {
total_filtered_reads = count_filter + count_other_filter
if (count_reads == 0) {
        pct_discarded_reads = 0
        pct_filtered_reads = 0
}else{
    pct_discarded_reads = count_discard* 100 / count_reads
    pct_filtered_reads = total_filtered_reads * 100 / count_reads
}

print "INDEX CORRECTION - SUMMARY STATS"
print "all reads:\t\t\t"count_reads
print "\tmatching reads:\t\t"count_match
print "\t\tfiltered:\t"count_filter
print "\t\tdiscarded:\t"count_discard
print "\t\tunchanged:\t"count_unchanged
print "\tother reads:\t\t"count_other
print "\t\tfiltered:\t"count_other_filter
print "\t\tunchanged:\t"count_other_unchanged

print "\nTotal filtered reads:\t\t"total_filtered_reads 
print "Percentage filtered reads:\t"pct_filtered_reads"%"
print "Percentage discarded reads:\t"pct_discarded_reads"%"
print "INDEX CORRECTION - SUMMARY STATS END"
}
' $index_file $map_file



#echo "map_file_index_correction.sh finished"