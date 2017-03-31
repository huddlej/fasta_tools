#!/bin/bash

# Convert TRF .dat output into a BED file.

if [[ "$#" -ne "1" ]]
then
    echo "Usage: $0 trf.dat > trf.bed"
    exit 1
fi

grep -P "^Sequence|^\d+" $1 \
    | sed 's/\://' \
    | awk 'OFS="\t" { if ($1 == "Sequence") { name=$2 } else { print name,$1,$2,$14,$2-$1 } }' \
    | sort -k 1,1 -k 2,2n
