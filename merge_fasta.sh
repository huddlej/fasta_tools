#!/bin/bash

# Usage: merge_fasta.sh newheader < one.fa two.fa three.fa > merged.fa

echo ">$1"

# Take two or more FASTA files from stdin, remove header lines and line returns,
# and split into lines of 60 characters long.
grep -vh "^>" | tr -d '\n' | fold -w 60
echo
