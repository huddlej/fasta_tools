#!/bin/bash

if [ $# -ne 2 ]
then
	echo "run under the context path directory. Prepare and format fasta files, prepare repeatmasker SGE file"
	echo "$0 accession_file outdir"
	echo
	echo "Please make sure the fasta files are named precisely after the header after this step"
	echo "The smooth working of the rest of the pipeline relys entirely on that"
	exit
fi

afile=$1
outdir=$2

# fasta download and format
mkdir -p ${outdir}
fasta_download  -m -f $afile -o ${outdir} &> fasta_download.log

cd ${outdir}
# clean header. Will keep the version number if there is any. But will be OK if no version is found
# in other words, it wil keep the name after gb that is only composed of letters, digits and dot
perl -i -pe "s/>gi\|\d+\|(gb|emb)\|([\w\.]+).+/>\2/" *

# rename file as header, in case version is not included at the begining
for fl in `find ./ -type f`
do
	nnm=`head -1 $fl | sed -re "s/>//"`
	mv $fl $nnm
  	echo -e "$fl\t$nnm" | sed 's/\.\///g'
done
cd ..