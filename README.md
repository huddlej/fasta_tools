# fasta_tools

Tools for working with FASTA files and related formats

Tool | Purpose
-----| -------
query_clones.py | Search NCBI for nucleotide accessions or sequences by clone names

## Query NCBI's nucleotide database by clone name

Search for clone accessions or sequences by clone name.

Usage:

```
usage: query_clones.py [-h] [--sequence] clone_names output

positional arguments:
  clone_names  name of a clone to search for in NCBI
  output       tab-delimited output of clone name/accession pairs OR sequence
               corresponding to clone names

optional arguments:
  -h, --help   show this help message and exit
  --sequence   download sequence corresponding to clone names instead of
               returning accessions
```
