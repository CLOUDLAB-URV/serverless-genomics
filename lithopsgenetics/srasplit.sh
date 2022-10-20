#!/bin/bash

name_seq=$1
num_spots=$2
db=$3
i=0
totalspots=0
div=0
mod=0
ini=0
end=0
chunk_length=0
bytes=0








let totalspots=$(awk -F "\"*,\"*" '{print $4}' info.csv | awk 'NR==2')
let bytes=$(awk -F "\"*,\"*" '{print $7}' info.csv | awk 'NR==2')

mod=$((totalspots % num_spots))


echo "Splitting files..."


PAIRED="$(awk -F "\"*,\"*" '{print $16}' info.csv | awk 'NR==2')"
VAR2="PAIRED"


echo "total spots:" $totalspots

while [ $end -lt $totalspots ]
do
#$((totalspots-ini)) -lt $num_spots
  if [ $end -eq $((totalspots-(num_spots+mod))) ]
  then
      let ini=$end
      let end=$totalspots
      echo $((ini)) $((end)) >> "${name_seq}.info"
  else
      if [  $i -eq 0 ]
        then
            let ini=$end
            let end=$end+$num_spots
            echo $((ini)) $((end))  >> "/tmp/${name_seq}.info"
        else
            let ini=$end
            let end=$end+$num_spots-1
            echo $((ini)) $((end))  >> "/tmp/${name_seq}.info"
        fi
      
        
    fi

    let i=$i+1
done

