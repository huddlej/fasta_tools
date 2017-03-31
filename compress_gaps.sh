#!/bin/bash
# Compress all runs of N's (gaps) into a single N each.

if [[ "$#" -ne "1" ]]
then
    echo "Usage: $0 input.fasta"
    exit 1
fi

INPUT=$1
INPUT_NAME=`basename ${INPUT}`

head -n 1 ${INPUT}
sed 1d ${INPUT} | tr -d '\n' | sed 's/N\+/N/g' | fold -w 80
