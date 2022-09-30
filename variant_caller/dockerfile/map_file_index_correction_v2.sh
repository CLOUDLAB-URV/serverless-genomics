#!/bin/bash
# script that takes as input:
# 1. filtered map file (with original line numbers, and best index field) as obtained from parse_gem_maxindex_minimapfile_stdin.sh
# 2. The corrected index file
# and removes any alignment that is below the threshold (max index value + tolerance) in terms of its index value.
# threshold is defined as "max index + tolerance" in the case of reads with no matches elsewhere, 
# or as "max corrected index + tolerance" for reads with matches in different chunks

echo "map_file_index_correction.sh starting"

# ARG 1: filtered map file 

map_file_default="/home/lumar/bioss/test_data/Hsapiens/map/SRR15068323_split_1_22__hg19_145450269split_00000001.se_filt_wline_no.map"
map_file=${1:-$map_file_default}

map_file_name="${map_file%.*}"

# ARG 2: corrected index file

index_file_default="/home/lumar/bioss/cloudbutton/scripts/07_aln_score_parser/merged_map_1_2_3_4_shared_only.txt"
index_file=$2:-$index_file_default}

# ARG 3: paired- v single-end sequencing
seq_type ${3: "single-end"}

# ARG 4: tolerance value
# specify how many strata to take below top index
tolerance=${4:-0}

# ARG 5: debug option
debug="${5:-"False"}"

echo -e "processing map file\n$map_file\nusing corrected index file \n$index_file\nwith a tolerance of **$tolerance** additional strata."

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
    # number of alignments to keep (above threshold)
    ok_aln = ""
    for (i = 1; i <= stratum_no; ++i) {
        if (strata[i]~"[1-9]") {
            if (i <= threshold) {
                ok_aln = ok_aln + strata[i]     #0:0:1:1:0  - define maximum index as $2 + tolerance
            }
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
    # match summary field
    summary=$6  
    # split match summary into two parts (before and after +)
    split(summary, summary_parts, /\+/)
    print_debug(debug,"match summary: "summary)
    print_debug(debug,"match summary part 1: "summary_parts[1])
    match_status = ""

    print_debug(debug,"index matching?: "$1"\t"s[1])
    # if line number matches, compare max corrected index (corrected index - s[2]) with max read index ($2) in chunk
    if ($1 == s[1]) { 
        print_debug(debug,"found match: "$1"\t"s[1])
        count_match++
        threshold = s[2] + TOLERANCE
        if ($2 <= threshold) {
            match_status="keep"
            # find number of alignments to keep based on threshold
            ok_aln = count_ok_aln(summary_parts[1], threshold)
            # find number of alignments listed in alignment field
            tot_aln = split($7, aligns, /,/)
            # if total alignments to keep < total number of alignments, then remove the unwanted alignments.
            if (tot_aln > ok_aln) {
                count_filter++
                filtered_aligns=filter_aln($7, ok_aln)
                #print t[$1]"\t"$1"\t"$2"\t"match_status"\t"summary"\t"ok_aln"\t"tot_aln"\t"$7
                #print t[$1]"\t"$1"\t"$2"\t"match_status"\t"summary"\t"ok_aln"\t"tot_aln"\t"filtered_aligns
                printf ("%s\t%s\t%s\t%s\t%s\n", $3, $4, $5, $6, filtered_aligns)  > NAME"_corrected.map"
            }
            else {
                count_unchanged++
                #print $0 > NAME"_corrected.map"
                printf ("%s\t%s\t%s\t%s\t%s\n", $3, $4, $5, $6, $7)  > NAME"_corrected.map"
            }
        }
        else {
            # discard read from chunk
            count_discard++
            match_status = "discard"
            print_debug(debug,"match below threshold: "t[$1]"\t"$1"\t"$2"\t"match_status"\t"summary)
            
            
        }
    }
    else {
        # if line number does not match, output as .map file omitting line number and index, and use threshold = "index + tolerance" 
        count_other++
        threshold = $2 + TOLERANCE
        ok_aln = count_ok_aln(summary_parts[1], threshold)
        # find number of alignments listed in alignment field
        tot_aln = split($7, aligns, /,/)
        # if total alignments to keep < total number of alignments, then remove the unwanted alignments.
        if (tot_aln > ok_aln) {
            count_other_filter++
            filtered_aligns=filter_aln($7, ok_aln)
            printf ("%s\t%s\t%s\t%s\t%s\n", $3, $4, $5, $6, filtered_aligns)  > NAME"_corrected.map"
        }
        else {
            count_other_unchanged++
            printf ("%s\t%s\t%s\t%s\t%s\n", $3, $4, $5, $6, $7)  > NAME"_corrected.map"
        }

    } 
}
}
END {
total_filtered_reads = count_filter + count_other_filter
pct_discarded_reads = count_discard* 100 / count_reads
pct_filtered_reads = total_filtered_reads * 100 / count_reads

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



echo "map_file_index_correction.sh finished"