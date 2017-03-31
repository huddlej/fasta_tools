#!/bin/bash
module load ucsc/20140617

if [[ "$#" -ne "2" ]]
then
    echo "Usage: $0 original.bb merged.bb"
    exit 1
fi

original=$1
original_name=$(basename $original)
intermediate_bed=/tmp/${original_name/.bb/.bed}
merged_bed=${intermediate_bed/.bed/.merged.bed}
output=$2

bigBedToBed ${original} ${intermediate_bed}
awk -f ~jlhudd/fasta_tools/merge_bed_by_copy.awk ${intermediate_bed} > ${merged_bed}
bedToBigBed ${merged_bed} /net/eichler/vol2/eee_shared/assemblies/hg19/chromInfo.txt ${output}

rm -f ${intermediate_bed} ${merged_bed}
