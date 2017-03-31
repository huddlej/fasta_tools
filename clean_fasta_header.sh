#!/bin/bash
#
# Replace FASTA header from NCBI with the accession.
#
# Usage: cat *.fasta | clean_fasta_header.sh > cleaned.fasta

#>gi|37537411|dbj|BS000144.1|

sed 's/>gi|[0-9]\+|\([a-zA-Z]\+\)|\([-_A-Za-z0-9\.]\+\).\+/>\2/'
