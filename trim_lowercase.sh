#!/bin/bash
# Trim all runs of lowercase letters from the beginning and the end of a sequence.

if [[ "$#" -ne "1" ]]
then
    echo "Usage: $0 input.fasta"
    exit 1
fi

INPUT=$1
INPUT_NAME=`basename ${INPUT}`

head -n 1 ${INPUT}
sed 1d ${INPUT} | tr -d '\n' | sed 's/[actgn]\+//g' | fold -w 80
echo
