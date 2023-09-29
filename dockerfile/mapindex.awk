BEGIN {
total_count=0
no_match_count=0
unlikely_count=0
filtered_count=0
}
{ 
    ## parse columns depending on single vs paired-end sequencing 
    if (seq_type == "paired-end" && NF<7) {
        #print_debug(debug, "alignment is not paired end - discarded")
        next
    }
    read_name=$1
    if (seq_type == "single-end") {
        read=$2
        qual=$3
        summary=$4
        aln_col=$5
    }
    else { # paired-end
        read1=$2
        read2=$3
        qual1=$4
        qual2=$5
        read=read1" "read2
        qual=qual1" "qual2
        summary=$6
        aln_col=$7
    }
    ## SELECT READS WITH A MATCH
    
    total_count++
	if (aln_col=="-") {
		no_match_count++
	}
    else {  

        #read_length=length(seq)

        # REMOVE UNLIKELY SECONDARY ALIGNMENTS
        # split summary on + to check for alignments vs secondary alignments
        # any reads with only secondary alignments are eliminated
        summary_info=split(summary, summary_parts, /\+/)
        if (summary_parts[1]~"[1-9]") {
            # RETRIEVE INDEX
            stratum_no = split(summary_parts[1], strata, /:/)
            idx = ""
            for (g = 1; g <= stratum_no; ++g) {
                #print strata[g]
                if (strata[g]~"[1-9]") {
                    idx = g
                    #print "index\t"g
                    break
                }
            }            
            ## print outputs
            # 1. map index file 
            printf ("%s\t%s\n", NR, idx) > NAME"_map.index.txt"
            # 2. filtered map file
            # with line number for correcting based on corrected index
            printf ("%s\t%s\t%s\n", NR, idx, $0)  > NAME"_filt_wline_no.map"
			filtered_count++

        }
		else {
		unlikely_count++
		}
    }
}
END {
print "total number of alignments lines: "total_count
print "alignment lines with no alignment: "no_match_count
print "unlikely alignments: "unlikely_count
print "filtered alignments: "filtered_count
count_sum = no_match_count + unlikely_count + filtered_count
print "sum of alignment counts: "count_sum
}
