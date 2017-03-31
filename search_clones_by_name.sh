#!/bin/bash

if [[ "$#" -ne "1" ]]
then
    echo "Usage: $0 clone_names.txt > clone_accessions_by_name.tab"
    exit 1
fi

clone_names=$1

while read clone_name
do
    accession=`python ~jlhudd/fasta_tools/search_clone_by_name.py ${clone_name}`;
    echo -e "${clone_name}\t${accession}";
done < ${clone_names}
