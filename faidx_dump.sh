#!/bin/bash

if [[ "$#" -lt "2" ]]
then
    echo "Usage: $0 reference region [region.fasta]"
    exit 1
fi

reference=$1
region=$2

if [[ "$#" -eq "3" ]]
then
    output=$3
else
    output="$(echo ${region} | sed 's/[-:]/_/g').fasta"
fi

samtools faidx ${reference} ${region} > ${output}
