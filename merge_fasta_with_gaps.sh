#!/bin/bash

if [[ "$#" -lt "1" ]]
then
    echo "Usage: $0 new_header one.fa two.fa three.fa > merged.fa"
    exit 1
fi

echo ">$1"
shift

for file in "$@"
do
    sed '/^>/d' ${file}
    cat ~jlhudd/projects/pacbio/smrtanalysis-2.2.0/gap.fasta
done | tr -d '\n' | sed 's/N\+$//' | fold -w 60

echo