#!/bin/bash

if [[ "$#" -lt "2" ]]
then
    echo "Usage: $0 <sequence path> <vector path>"
    exit 1
fi

sequence_path=$1
vector_path=$2
tmp_dir=$3

sequence_name=`basename ${sequence_path}`

if [[ -z "${tmp_dir}" ]]
then
    tmp_dir="/tmp/${sequence_name}"
fi

mkdir -p ${tmp_dir}

bl2seq -p blastn -m T -F F -i ${vector_path} -j ${sequence_path} -o ${tmp_dir}/vector.bo -D 1

# Get vector's coordinates in the contig, sort, and merge. Take the largest
# region after merging as the vector location.
grep -v "^#" ${tmp_dir}/vector.bo \
    | cut -f 2,9,10 \
    | awk 'OFS="\t" { if ($2 > $3) { print $1,$3,$2 } else { print $0 } }' \
    | sort -k 1,1 -k 2,2n \
    | mergeBed -i stdin -d 1 \
    | awk 'OFS="\t" { print $0,$3-$2 }' \
    | sort -k 4,4rn \
    | head -n 1 > ${tmp_dir}/vector.bed

lines=`wc -l ${tmp_dir}/vector.bed | awk '{ print $1 }'`

if [[ "${lines}" -gt "0" ]]
then
    fasta_length -o ${tmp_dir}/fasta_length.log ${sequence_path} &> /dev/null
    sed '1d;$d' ${tmp_dir}/fasta_length.log | cut -f 3-4 > ${tmp_dir}/lengths.tab
    complementBed -i ${tmp_dir}/vector.bed -g ${tmp_dir}/lengths.tab > ${tmp_dir}/contigs.bed
    bedtools getfasta -fi ${sequence_path} -bed ${tmp_dir}/contigs.bed -fo ${tmp_dir}/${sequence_name}.tab -tab
    bedtools getfasta -fi ${sequence_path} -bed ${tmp_dir}/contigs.bed -fo ${tmp_dir}/${sequence_name}.fasta

    sed -i 's/[:]/\t/g;s/\t\([0-9]\+\)-\([0-9]\+\)/\t\1\t\2/' ${tmp_dir}/${sequence_name}.tab
    sort -k 1,1 -k 2,2rn ${tmp_dir}/${sequence_name}.tab | cut -f 4 \
        | ~jlhudd/fasta_tools/merge_fasta.sh ${sequence_name} > `pwd`/${sequence_name}.no_vector.fasta
    cp ${tmp_dir}/${sequence_name}.fasta `pwd`/${sequence_name}.pieces.fasta
else
    cp ${sequence_path} `pwd`/${sequence_name}.no_vector.fasta
fi

#rm -rf ${tmp_dir}
