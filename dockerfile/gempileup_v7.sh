#!/bin/bash
# script that converts a .map alignment file to mpileup
# MAP FIELDS:
# $1 = READ NAME
# $2 = SEQUENCE(S)
# $3 = QUALITY SCORES
# $4 = MATCH SUMMARY
# $5 = ALIGNMENTS

# legend:
# indel+ = deletion in genome
# indel- = deletion in read. 
# base_count = counting bases in read as parsing gigar
# read_start = start position of read match on chromosome, indicating the first base on the forward strand
# abs_bpos = chromosome position of each base called in the mpileup output
# rel_bpos = relative position of each base called within the trimmed read

# test commandline (2)
# debug command-line
# bash gempileup_v7.sh SRR15068323_split_1_22__hg19_145450269split_00000001.se.map /home/lumar/bioss/test_data/Hsapiens/fasta/hg19_145450269split_00000001.fasta True > SRR15068323_split_1_22__hg19_145450269split_00000001.se.map.nomerge.debug.mpileup
# commandline for mpileup output
# ./cloudbutton/scripts/10_gempileup/SRR15068323_split_1_22__hg19_145450269split_00000001.se.map
# bash gempileup_v7.sh SRR15068323_split_1_22__hg19_145450269split_00000001.se.map /home/lumar/bioss/test_data/Hsapiens/fasta/hg19_145450269split_00000001.fasta   | bash  gempileup_merge.sh > SRR15068323_split_1_22__hg19_145450269split_00000001.se.map.mpileup
# awk  '{A[$2]++}END{for(i in A) {if(A[i]>1) {print i,A[i]}}}' SRR15068323_split_1_22__hg19_145450269split_00000001.se.map.mpileup | sort -k 1,1n -k 2,2n > duplicates20.txt
#awk 'FNR>=5 && FNR<=9' /home/lumar/bioss/test_data/Hsapiens/map/SRR15068323_split_1_22__hg19_145450269split_00000001.se_filt_wline_no_corrected.map
#file=$1
# filtered and corrected map file 
#file1="/home/lumar/bioss/test_data/Hsapiens/map/SRR15068323_split_1_22__hg19_145450269split_00000001.se_filt_wline_no_corrected.map"
#echo "processing "$file1
#file1_name="${file1%.*}"

# <com> usage example
#cat MycobacteriumTuberculosis.CompleteGenomes.fasta | fasta-tabular | awk -F '\t' 'BEGIN{rc="rev|com"; com="com"} {print $2 |& com; com |& getline line; print ">"$1"\n"line}' > TEST

map_file="${1:-/dev/stdin}"
ref_genome="${2:-/dev/stdin}"
debug="${3:-"False"}"
#echo "gempileup_v7.sh starting"
awk -v NAME="$map_file" -v ref_genome="$ref_genome" -v debug="$debug" '
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
function print_newline_no_debug(debug)
{
    if (debug !="True") {
        printf("\n")
    }
}
function join(array, start, end, sep, result, i)  #https://www.gnu.org/software/gawk/manual/html_node/Join-Function.html
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
function complement(arg)
{
    print arg |& com; com |& getline compl_string
    return compl_string
}
#function remove_clip(gig_el,base_count,n,start,seq)
function remove_clip(seq)
{
    clip_n = gensub(/\(|\)/, "","g",gig_el[n]) 
    print_debug(debug,"clip no. "gig_el[n]" length: "clip_n)
    # reset the start variable in case there is a soft clip at the beginning of the read
    if (n==1) { #remove sequence from beginning of read
        start=start-1
        new_seq_start = clip_n + 1
        print_debug(debug,"incoming sequence: "seq)
        print_debug(debug,"new sequence start: "new_seq_start)
        seq = substr(seq, new_seq_start)
        print_debug(debug,"new sequence: "seq)
    }
    return(seq)
}

#function process_mutation_indel(read_info,gigar_index,strand,seq,qual_string,gig_elem,base_count,start,base_out,abs_bpos,read_start,ins,del,ins_prev,del_prev,ins_cumul,del_cumul,ins_cumul_prev,del_cumul_prev,chr,read_len) 
function process_mutation_indel(seq,gig_elem) 
{
    base_count = base_count + 1
    
    # 1. process mutation
    if (strand=="+") {
        abs_bpos = read_start + base_count -1 + del_cumul
        ref_bp = substr(gig_elem,1,1) 
        base_out_mut = substr(seq,base_count+ins_cumul,1)
        if (start==1) {
            base_out_mut = "^"substr(seq,base_count+ins_cumul,1)
        }
        qual = substr(qual_string,base_count+ins_cumul,1)
    }
    else {
        abs_bpos = read_start + read_len - base_count - del_cumul
        ref_bp = complement(substr(gig_elem,1,1)) 
        base_out_mut = tolower(substr(seq,base_count+ins_cumul,1))
        if (start==1) {
            base_out_mut = "^"tolower(substr(seq,base_count+ins_cumul,1))
        }
        qual = substr(qual_string,base_count+ins_cumul,1)
    }

    # 2. add indel information to mutation line, if indel present
    if (gig_elem~/>/) { # indel present in element
        #print_debug(debug, "pre-number parsing seq: "seq)
        # separate number of matches to reference from indel information (i.e. split 40>>5+ to 40 and >5+)
        split(gig_elem, n_and_indel, />/)
        match_elem = n_and_indel[1]
        indel_gigar = n_and_indel[2]
        indel_type = substr(indel_gigar,length(indel_gigar),1)
        # extract indel length
        #indel_length = gensub(/>|indel_type/, "","g",indel_gigar) 
        indel_length = int(gensub(/[+-]/, "","g",indel_gigar)) 
        print_debug(debug, indel_type"indel ~ length: "indel_length)
        
        # assign indel length to insertion / deletion variable
        if (indel_type=="+") { #deletion in genome
            #print "checking before - del: "del"; del_cumul: "del_cumul
            del = indel_length
            del_cumul = del_cumul_prev + del
            #print "checking after - del: "del"; del_cumul: "del_cumul
            del_seq=toupper(substr(seq_array[chr], abs_bpos+1, indel_length))
            #base_out_mut=base_out_mut"-"indel_length""del_seq
            }
        else if (indel_type=="-") { # insertion in read
            #print "checking before - ins: "ins"; ins_cumul: "ins_cumul"; ins_cumul_prev: "ins_cumul_prev
            ins = indel_length
            ins_cumul = ins_cumul_prev + ins
            #print "checking after - ins: "ins"; ins_cumul: "ins_cumul"; ins_cumul_prev: "ins_cumul_prev
            ins_seq=substr(seq,rel_bpos+1,indel_length)
            #base_out_mut=base_out_mut"+"indel_length""ins_seq
        }
    }

    # 3. print mutation line
    printf_debug(debug,"%s\t%s%s\t", read_info, gigar_index, strand) 
    printf ("%s\t%s\t%s\t%s\t%s\t%s\t", chr, abs_bpos, ref_bp, "1", base_out_mut, qual)
    print_newline_no_debug(debug)
    print_debug(debug, "read_start: "read_start"; rel_bpos: "rel_bpos"; base_count: "base_count"; m: "m"; del: "del"; ins: "ins"; del_prev: "del_prev"; ins_prev: "ins_prev"; del_cumul: "del_cumul"; ins_cumul: "ins_cumul"; del_cumul_prev: "del_cumul_prev"; ins_cumul_prev: "ins_cumul_prev)

    # 4. update variables
    ins_cumul_prev = ins_cumul
    del_cumul_prev = del_cumul
    #print "checking after end of element loop - ins: "ins"; ins_cumul: "ins_cumul"; ins_cumul_prev: "ins_cumul_prev
    del_prev = del
    ins_prev = ins
    
    # 5. print deleted positions in genome on separate lines
    if (indel_type=="+") { # deleted positions in genome
        for (nm = 1; nm <= del; nm++) {
            if (strand == "+") {
                abs_bpos++
            }
            else if (strand == "-") {
                abs_bpos--
            }
            printf_debug(debug,"%s\t%s%s\t", read_info, gigar_index, strand) 
            base_from_ref = toupper(substr(seq_array[chr], abs_bpos, 1))
            printf ("%s\t%s\t%s\t%s\t%s\t%s\t", chr, abs_bpos, base_from_ref, "1", "*", "*")
            print_newline_no_debug(debug)
            print_debug(debug, "read_start: "read_start"; rel_bpos: "rel_bpos"; base_count: "base_count"; m: "m"; del: "del"; ins: "ins"; del_prev: "del_prev"; ins_prev: "ins_prev"; del_cumul: "del_cumul"; ins_cumul: "ins_cumul"; del_cumul_prev: "del_cumul_prev"; ins_cumul_prev: "ins_cumul_prev)
        }
    }

    # 6. reset variables
    ins=0
    del=0 

    return(base_count"_"ins_prev"_"del_prev"_"ins_cumul"_"del_cumul"_"ins_cumul_prev"_"del_cumul_prev)
}

#function process_number_indel(read_info,gigar_index,strand,seq,qual_string,gig_elem,base_count,start,base_out,abs_bpos,read_start,ins,del,ins_prev,del_prev,ins_cumul,del_cumul,ins_cumul_prev,del_cumul_prev,chr,read_len) 
function process_number_indel(seq,gig_elem) 

{
    if (gig_elem!~/>/) {
        for (m = 1; m <= gig_elem; m++) {
            rel_bpos = base_count + m + ins_cumul_prev
            if (strand=="+") {
                abs_bpos = read_start + base_count + m - 1 + del_cumul
            }
            else {
                abs_bpos = read_start + read_len - base_count - m - del_cumul
            }
            ref_bp = substr(seq,rel_bpos,1)
            qual = substr(qual_string,rel_bpos,1)
            if (start==1 && (m-base_count)==1) {
                base_out="^"base_out
                qual="I"qual
            }
            printf_debug(debug,"%s\t%s%s\t", read_info, gigar_index, strand) 
            printf ("%s\t%s\t%s\t%s\t%s\t%s\t", chr, abs_bpos, ref_bp, "1", base_out, qual)
            print_newline_no_debug(debug)
            print_debug(debug, "read_start: "read_start"; rel_bpos: "rel_bpos"; base_count: "base_count"; m: "m"; del: "del"; ins: "ins"; del_prev: "del_prev"; ins_prev: "ins_prev"; del_cumul: "del_cumul"; ins_cumul: "ins_cumul"; del_cumul_prev: "del_cumul_prev"; ins_cumul_prev: "ins_cumul_prev)
            
            # reset values after printing
            if (strand=="+") {
                base_out="."
            }
            else {
                base_out=","
            }
        }
        base_count = base_count + gig_elem
        # reset variables
        m = 0

    }
    else {  # indel present in element
        #print_debug(debug, "pre-number parsing seq: "seq)
        # separate number of matches to reference from indel information (i.e. split 40>>5+ to 40 and >5+)
        split(gig_elem, n_and_indel, />/)
        match_elem = n_and_indel[1]
        indel_gigar = n_and_indel[2]
        indel_type = substr(indel_gigar,length(indel_gigar),1)
        # extract indel length
        indel_length = gensub(/>|indel_type/, "","g",indel_gigar) 
        indel_length = int(gensub(/[+-]/, "","g",indel_length)) 
        print_debug(debug, indel_type"indel ~ length: "indel_length)
        
        # assign indel length to insertion / deletion variable
        if (indel_type=="+") { #deletion in genome
            #print "checking before - del: "del"; del_cumul: "del_cumul
            del = indel_length
            del_cumul = del_cumul_prev + del
            #print "checking after - del: "del"; del_cumul: "del_cumul
            }
        else {  # insertion in read
            #print "checking before - ins: "ins"; ins_cumul: "ins_cumul"; ins_cumul_prev: "ins_cumul_prev
            ins = indel_length
            ins_cumul = ins_cumul_prev + ins
            #print "checking after - ins: "ins"; ins_cumul: "ins_cumul"; ins_cumul_prev: "ins_cumul_prev
        }


        # parse matches and include indel information in last match in element
        for (m = 1; m <= match_elem; m++) {

            rel_bpos = base_count + m + ins_cumul_prev
            ref_bp = substr(seq,rel_bpos,1)
            qual = substr(qual_string,rel_bpos,1)
            if (strand=="+") {
                abs_bpos = read_start + base_count + m - 1 + del_cumul_prev
            }
            else {
                abs_bpos = read_start + read_len - base_count - m - del_cumul_prev
            }
            
            # format first match in element
            if (start==1 && (m-base_count)==1) {
                base_out="^"base_out
                qual="I"qual
            }

            # format last match in element
            if (m == match_elem) { # print indel information in last match in element
                # indel printing removed temporarily
                if (indel_type=="+") {  #deletion in genome
                    #mpileup-formatted output ("If there is a deletion after this read base, text matching “-[0-9]+ [ACGTNacgtn]+”: a “-” character followed by the deleted reference bases represented similarly." -mpileup specifications)
                    if (strand=="+") { 
                        indel_seq=substr(seq_array[chr], abs_bpos+1, indel_length)
                    }
                    else {
                        indel_seq=substr(seq_array[chr], abs_bpos-indel_length, indel_length)
                    }
                    base_out=base_out"-"indel_length
                }
                else if (indel_type=="-") { # insertion in read
                    #mpileup-formatted output (uses "+" for insertion in read, while .map file uses "-")
                    indel_seq=substr(seq,rel_bpos+1,indel_length)
                    base_out=base_out"+"indel_length""indel_seq
                }
            }

            # print mpileup line
            printf_debug(debug,"%s\t%s%s\t", read_info, gigar_index, strand) 
            printf ("%s\t%s\t%s\t%s\t%s\t%s\t", chr, abs_bpos, ref_bp, "1", base_out, qual)
            print_newline_no_debug(debug)
            print_debug(debug, "read_start: "read_start"; rel_bpos: "rel_bpos"; base_count: "base_count"; m: "m"; del: "del"; ins: "ins"; del_prev: "del_prev"; ins_prev: "ins_prev"; del_cumul: "del_cumul"; ins_cumul: "ins_cumul"; del_cumul_prev: "del_cumul_prev"; ins_cumul_prev: "ins_cumul_prev)
            # remove "^" after assigning it to first base of read.
            if (base_out~/^\^/) {
                base_out=gensub(/\^/, "","g",base_out) 
            }

        }
        # update base count
        base_count = base_count + match_elem 
        # reset match iteration
        m = 0
        # reset base_out variable after writing gigar string containing indel
        if (strand=="+") { 
            base_out="."
        }
        else {
            base_out=","
        }
        ins_cumul_prev = ins_cumul
        del_cumul_prev = del_cumul
        #print "checking after end of element loop - ins: "ins"; ins_cumul: "ins_cumul"; ins_cumul_prev: "ins_cumul_prev
        del_prev = del
        ins_prev = ins
        

        # print deleted positions in genome
        if (indel_type=="+") { # deleted positions in genome
            for (nm = 1; nm <= del; nm++) {
                if (strand == "+") {
                    abs_bpos++
                }
                else if (strand == "-") {
                    abs_bpos--
                }
                printf_debug(debug,"%s\t%s%s\t", read_info, gigar_index, strand) 
                base_from_ref = toupper(substr(seq_array[chr], abs_bpos, 1))
                printf ("%s\t%s\t%s\t%s\t%s\t%s\t", chr, abs_bpos, base_from_ref, "1", "*", "*")
                print_newline_no_debug(debug)
                print_debug(debug, "read_start: "read_start"; rel_bpos: "rel_bpos"; base_count: "base_count"; m: "m"; del: "del"; ins: "ins"; del_prev: "del_prev"; ins_prev: "ins_prev"; del_cumul: "del_cumul"; ins_cumul: "ins_cumul"; del_cumul_prev: "del_cumul_prev"; ins_cumul_prev: "ins_cumul_prev)
            }
        }
        ins=0
        del=0 
    }
    #print_debug(debug, "end of number/indel function - base count: "base_count"; ins: "ins"; del: "del)
    return(base_count"_"ins_prev"_"del_prev"_"ins_cumul"_"del_cumul"_"ins_cumul_prev"_"del_cumul_prev)
}

BEGIN {
    # com - executable for generating complementary DNA bases
    com = "com"
    count_reads=0 # all reads
    num=0
    while (getline < ref_genome) { 
        if ($0~"^>"){
            num++
            #remove trailing whitespace
            fasta_header=$0
            fasta_header=gensub(/>/, "","g",fasta_header)
            fasta_header=gensub(/[ \t]+$/, "","g",fasta_header)
            #print "header no.: "num": "$1
            
            # remove > sign
            #print "fasta header with > removed: "fasta_header
        }
        else{
            #print "sequence detected, associated with fasta header: "fasta_header
            seq_array[fasta_header]=$0
            #print "seq "num":"$0
            #print "seq from array: " seq_array[fasta_header]
            #print "substring test header 1: "substr(fasta_header, 1, 2)
            #print "substring test seq 1: "substr(seq_array[fasta_header], 12000, 10)
        }
    }
    
    #for (i in seq_array) print i"\t"seq_array[i]
    #print "substring test 2 header: "substr(fasta_header, 1, 2)
    #print "substring test 2 seq: "substr(seq_array[fasta_header], 12000, 10)
}
{ 
count_reads++
#if (count_reads >=1 && count_reads <=23000 && $5!="-") {
if (count_reads >= 1 && $5!="-") {
#if ($0~"SRR15068323.9174_") {
    seq = $2
    read_len = length($2)
    print_debug(debug,"read length first: "read_len)
    qual_string = $3

    # if alignments are on minus strand, create complement (no need to reverse qual_string)
    compl_seq=""
    if ($5~":-:") { # check for minus strand in alignment field
        compl_seq = complement(seq)
    }
    
    # SPLIT ALIGNMENTS FIELD INTO INDIVIDUAL ALIGNMENTS
    align_n = split($5, aligns, /,/)
    print_debug(debug, "####\n##read "count_reads"\t"$5)
    #print_debug(debug,"number of alignments: "align_n)
    for (f = 1; f <= align_n; f++) {

        # split each alignment into columns (chr, strand, read_start, gigar)
        split(aligns[f], align_cols, /:/)
        chr=align_cols[1]
        strand=align_cols[2]
        read_start=align_cols[3]
        gigar=align_cols[4]
        # reset insertion and deletion values
        ins=0
        del=0 
        ins_prev=0
        del_prev=0
        ins_cumul=0
        del_cumul=0
        ins_cumul_prev=0
        del_cumul_prev=0

        # filter alignments containing double indels (i.e. >5+>7-)
        # go to next alignment
        if (gigar~/>[0-9]+[\+-]>[0-9]+[\+-]/) {
            continue
        }
        

        # MODIFY GIGAR TO ALLOW SPLITTING 
        gigar_mod=gensub(/([A-Za-z])/, "~\\1","g",gigar) 
        gigar_mod=gensub(/([0-9]+~)/, "~\\1","g",gigar_mod) 
        gigar_mod=gensub(/([0-9]+$)/, "~\\1","g",gigar_mod) 
        gigar_mod=gensub(/([0-9]+)(>)/, "~\\1\\2","g",gigar_mod) 
        gigar_mod=gensub(/([0-9]+)(\()/, "~\\1~\\2","g",gigar_mod) 
        gigar_mod=gensub(/^~/, "","g",gigar_mod) 
        read_info=NR":"f"/"align_n
        print_debug(debug,"#aln: "NR":"f"/"align_n"\t"strand"\t"gigar_mod"\t"$2)
        print_debug(debug, "gigar: "aligns[f])

        # SPLIT MODIFIED GIGAR
        split(gigar_mod, gig_el, "~", gig_char)

        #set read bp strand info based on alignment strand info
        # base output: nucleotide mutation, "." or ","
        base_out = ""
        
        # PLUS STRAND
        if (strand=="+") {  # alignment on + strand
            base_out="."
            base_count = 0
            start = 0  # variable to specify start of read
            seq = $2
            # ITERATE THROUGH GIGAR ELEMENTS (+ STRAND)
            for (n = 1; n <= length(gig_el); n++) {
                start++ 
                gigar_index=n"/"length(gig_el)"gig"
                print_debug(debug, "\t"n"/"length(gig_el)"\t"gig_el[n]"\t"gigar_mod)
                # PARSE DIFFERENT OPTIONS
                # 1. GIGAR ELEMENT IS A SOFT CLIP
                if (gig_el[n]~/\(/) {
                    seq = remove_clip(seq)
                    qual_string = remove_clip(qual_string)
                }
                # 2. GIGAR ELEMENT IS A NUMBER (AND POTENTIALLY CONTAINS INDEL)
                else if (gig_el[n]~/^[0-9]+/) {
                    base_count__ins__del = process_number_indel(seq,gig_el[n]) 
                    split(base_count__ins__del, base_count__ins__del_el, "_")
                    base_count = base_count__ins__del_el[1]
                    ins_prev = base_count__ins__del_el[2]
                    del_prev = base_count__ins__del_el[3]
                    ins_cumul = base_count__ins__del_el[4]
                    del_cumul = base_count__ins__del_el[5]
                    ins_cumul_prev = base_count__ins__del_el[6]
                    del_cumul_prev = base_count__ins__del_el[7]
                }
                # 3. GIGAR ELEMENT IS A MUTATION
                else if (gig_el[n]~/[A-Za-z]/) {
                    base_count__ins__del = process_mutation_indel(seq,gig_el[n]) 
                    # split function return value
                    split(base_count__ins__del, base_count__ins__del_el, "_")
                    # update count and position
                    base_count = base_count__ins__del_el[1]
                    ins_prev = base_count__ins__del_el[2]
                    del_prev = base_count__ins__del_el[3]
                    ins_cumul = base_count__ins__del_el[4]
                    del_cumul = base_count__ins__del_el[5]
                    ins_cumul_prev = base_count__ins__del_el[6]
                    del_cumul_prev = base_count__ins__del_el[7]
                }
            }
        }
        # MINUS STRAND
        else if (strand=="-") { # alignment on - strand
            base_out=","
            # go through each gigar element
            base_count = 0
            print_debug(debug, "complement seq: "compl_seq)
            print_debug(debug, "read length: "read_len)

            # modify read length
            tot_clip_n=0
            indel_plus=0
            indel_minus=0
            for (n = 1; n <= length(gig_el); n++) {
                if (gig_el[n]~/\(/) {
                    clip_n = gensub(/\(|\)/, "","g",gig_el[n])
                    tot_clip_n=tot_clip_n + clip_n 
                }
                else if (gig_el[n]~/>/) {
                    indel=0
                    #print gig_el[n]
                    split(gig_el[n], n_and_indel, />/)
                    #print "indel unedited: "n_and_indel[2]
                    # extract indel information
                    if (n_and_indel[2]~/-/) { # skip positions in read
                        indel = gensub(/[>-]/, "","g",n_and_indel[2]) 
                        print_debug(debug, "ins:-"indel"-")
                        print_debug(debug, "indel_minus before addition: "indel_minus)
                        print_debug(debug, "type of indel_minus: "typeof(indel_minus))
                        print_debug(debug, "type of indel: "typeof(indel))
                        print_debug(debug, "type of indel (int): "typeof(int(indel)))
                        indel_minus = indel_minus + int(indel)
                        print_debug(debug, "type of indel_minus: "typeof(indel_minus))
                        print_debug(debug, "indel_minus after addition: "indel_minus)
                    }
                    else if (n_and_indel[2]~/+/) {  # skip positions in genome (deletion in genome)
                        indel = gensub(/[>+]/, "","g",n_and_indel[2]) 
                        print_debug(debug, "del:-"indel"-")
                        print_debug(debug, "indel_plus before addition: "indel_plus)
                        print_debug(debug, "type of indel_plus: "typeof(indel_plus))
                        print_debug(debug, "type of indel: "typeof(indel))
                        print_debug(debug, "type of indel (int): "typeof(int(indel)))
                        indel_plus = indel_plus + int(indel)
                        print_debug(debug, "type of indel_plus: "typeof(indel_plus))
                        print_debug(debug, "indel_plus after addition: "indel_plus)
                    }
                }
            }
            base_n_adjust = indel_plus - tot_clip_n - indel_minus
            read_len = read_len + base_n_adjust 
            print_debug(debug, "total soft clip: "tot_clip_n)
            print_debug(debug, "total indel -: "indel_minus)
            print_debug(debug, "total indel +: "indel_plus)
            print_debug(debug, "total trimmed bases to add (+) / remove (-): "base_n_adjust)
            print_debug(debug, "modified read length: "read_len)

            # ITERATE THROUGH GIGAR ELEMENTS (- STRAND)
            for (n = 1; n <= length(gig_el); n++) {
                start++
                print_debug(debug, "\t"n"/"length(gig_el)"\t"gig_el[n]"\t"gigar_mod)
                # PARSE DIFFERENT OPTIONS
                # 1. GIGAR ELEMENT IS A SOFT CLIP
                if (gig_el[n]~/\(/) {
                    compl_seq = remove_clip(compl_seq)
                    qual_string = remove_clip(qual_string)
                }
                # 2. GIGAR ELEMENT IS A NUMBER (AND POTENTIALLY CONTAINS INDEL)
                else if (gig_el[n]~/^[0-9]+/) {
                    # process number (and indel) with function
                    base_count__ins__del = process_number_indel(compl_seq,gig_el[n]) 
                    # split function return value
                    split(base_count__ins__del, base_count__ins__del_el, "_")
                    # update count and position
                    base_count = base_count__ins__del_el[1]
                    ins_prev = base_count__ins__del_el[2]
                    del_prev = base_count__ins__del_el[3]
                    ins_cumul = base_count__ins__del_el[4]
                    del_cumul = base_count__ins__del_el[5]
                    ins_cumul_prev = base_count__ins__del_el[6]
                    del_cumul_prev = base_count__ins__del_el[7]
                }
                # 3. GIGAR ELEMENT IS A MUTATION
                else if (gig_el[n]~/[A-Za-z]/) {
                    base_count__ins__del = process_mutation_indel(compl_seq,gig_el[n]) 
                    # split function return value
                    split(base_count__ins__del, base_count__ins__del_el, "_")
                    # update count and position
                    base_count = base_count__ins__del_el[1]
                    ins_prev = base_count__ins__del_el[2]
                    del_prev = base_count__ins__del_el[3]
                    ins_cumul = base_count__ins__del_el[4]
                    del_cumul = base_count__ins__del_el[5]
                    ins_cumul_prev = base_count__ins__del_el[6]
                    del_cumul_prev = base_count__ins__del_el[7]
                }
            }
        }
    read_len = length($2)
    seq = $2
    compl_seq = complement(seq)
    }
} 
}
END {
}
' $map_file

#echo "gempileup_v7.sh completed"
