#!/bin/bash
# merge gempileup.sh output by position
#!/bin/bash
file="${1:-/dev/stdin}"
#echo "gempileup_merge.sh starting"
sort --parallel 3 -T . -k 1,1 -k 2,2n | 
awk -F '\t' '
function print_current()
    {
        if (old!="") print old"\t"len"\t"syms"\t"quals
    } 
{curr=$1"\t"$2"\t"$3; 
    if (curr!=old) {
        if (old!="") 
            print_current(); 
            len=0; 
            syms=""; 
            quals=""; 
            old=curr
        } 
        len+=$4; 
        syms=syms $5; 
        quals=quals $6
    } 
END{print_current()}
' $file
#echo "gempileup_merge.sh completed"


