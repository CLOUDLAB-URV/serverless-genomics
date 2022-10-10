#!/bin/bash
# script that provides the best match for each read.

# Legend:
# 1: line number
# 2: index

file1=$1
file2=$2 

awk '
{ 
if (NR==FNR) {
        t[$1]=$0
        coln_1 = NF
}
else {
    coln_2 = NF
    split(t[$1],s,"\t");
        #print $1"\t"s[1]
        # if line number matches select the best index (lowest number)
        best=""
        if ($1 == s[1]) { 
            if ($2 <= s[2]) {
                best = $2
            }    
            else {
                best = s[2]
            }
        flag="2F";
        printf ("%s\t%s\t%s\t%s\n", s[1],best, flag, "") 
        # delete array key (file 1)
        delete t[$1] 
        }
        else {
            # keep original flag (might be either 1F or 2F depending on previous match)
            choice="2U"
            if (coln_2==2) {
                printf ("%s\t%s\t%s\t%s\n", $1, $2, "1F", choice) 
            }
            else {
                if (coln_2==3||coln_2==4) {
                    printf ("%s\t%s\t%s\t%s\n", $1, $2, $3, choice) 
                }
            }
        } 

    }
}

END{
    # keys left in a were not in file2 (and thus were not deleted), so print their lines unchanged
    for (k in t){
        if (length(t[k])>0) {
            choice="1U"
            if (coln_1==2) {
                printf ("%s\t%s\t%s\n", t[k], "1F", choice) 
            }
            else {
                if (coln_1==3) {
                    printf ("%s\t%s\t%s\n", t[k], choice) 
                }
                else {
                    if (coln_1==4) {
                        printf (t[k]"_1U\n") 
                    }
                }
            }
        }
    }
} 
' $file1 $file2

