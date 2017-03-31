CONTIG="$1"

fasta_findNs.pl -i "$CONTIG"
cut -f 1-3 Npositions.tbl | sed 1d | awk 'OFS="\t"{ print $1,$2-1,$3 }' > gaps.bed

fasta_length -o length.log "$CONTIG"
cut -f 2,4 length.log | sed '1d;$d' > tmp && mv -f tmp length.log

complementBed -i gaps.bed -g length.log > regions.bed
fastaFromBed -fi "$CONTIG" -bed regions.bed -fo "$CONTIG.split"
sed -i 's/:/_/g;s/-/_/g' "$CONTIG.split"
~jlhudd/fasta_tools/split_fasta_file.py -s "" -r ">(.+)" "$CONTIG.split"

rm -f regions.bed length.log gaps.bed Npositions.tbl *.split *.fai