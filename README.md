# fasta_tools

Tools for working with FASTA files and related formats

Tool | Purpose
-----| -------
query_clones.py | Search NCBI for nucleotide accessions or sequences by clone names
fragment_fasta_records.py | Fragment FASTA records into sequences of a fixed size with optional sliding window

## Query NCBI's nucleotide database by clone name

Search for clone accessions or sequences by clone name.

Usage:

```
usage: query_clones.py [-h] [--sequence] clone_names output

positional arguments:
  clone_names  text file containing names of clones to search for in NCBI (one
               name per line)
  output       tab-delimited output of clone name/accession pairs OR FASTA
               sequence corresponding to clone names

optional arguments:
  -h, --help   show this help message and exit
  --sequence   download sequence corresponding to clone names instead of
               returning accessions
```

## Fragment FASTA sequence by a fixed size

Usage :

```
usage: fragment_fasta_records.py [-h] [--slide SLIDE] input output window

positional arguments:
  input          FASTA sequences to fragment
  output         fragmented FASTA sequences
  window         length to fragment each FASTA sequence to

optional arguments:
  -h, --help     show this help message and exit
  --slide SLIDE  length to fragment each FASTA sequence to
```
