#!/bin/env python
"""
Given a clone name, report any matching accessions in NCBI's nucleotide
database.
"""
import argparse
from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import csv
import sys
import urllib2
Entrez.email = "anonymous@anonymous.org"

# Set maximum records to return from NCBI's Entrez e-search to the maximum
# allowed by the service.
MAX_RECORDS_TO_RETURN = 50000


def get_ids_by_accession(accessions):
    """
    Search NCBI's nucleotide database for the given accessions and return a list
    of nucleotide database ids.

    >>> get_ids_by_accession(["AC254980.1"])
    ['595582430']
    >>> get_ids_by_accession(["fakeclone"])
    []
    """
    # By default, expect no results.
    results = []

    # Process clones in batches to avoid overly long URL requests.
    total = len(accessions)
    step = 100
    step_range = range(0, total, step)

    for start in step_range:
        end = min(start + step, total)
        step_accessions = accessions[start:end]

        # Search NCBI nucleotide database by accession.
        handle = Entrez.esearch(db="nucleotide", term=" OR ".join(step_accessions), retmax=MAX_RECORDS_TO_RETURN)
        search_record = Entrez.read(handle)

        if int(search_record["Count"]) > 0:
            results.extend(search_record["IdList"])

    # Get all distinct ids.
    results = list(set(results))

    return results


def download_sequences(nucleotide_ids):
    """
    Download sequences from NCBI's nucleotide database for the given nucleotide
    database ids and return a list of BioPython SeqRecords.

    >>> download_sequences(["146189515"])
    [SeqRecord(seq=Seq('GCGCCCAATACGCAAACCGCCTCTCCCCGCGCGTTGGCCGATTCATTAATGCAG...AGA', Alphabet()), id='AM697672.1', name='<unknown name>', description='Cloning vector pZCAA', dbxrefs=[])]
    >>> download_sequences(["fakeid"])
    []
    """
    output_records = []

    try:
        handle = Entrez.efetch(db="nuccore", id=",".join(nucleotide_ids), rettype="fasta", retmode="xml", retmax=MAX_RECORDS_TO_RETURN)
        nuccore_records = Entrez.read(handle)

        for nuccore_record in nuccore_records:
            output_records.append(SeqRecord(Seq(nuccore_record["TSeq_sequence"]), id=nuccore_record["TSeq_accver"], description=nuccore_record["TSeq_defline"]))
    except urllib2.URLError, e:
        sys.stderr.write("Error downloading sequences: %s\n" % e)

    return output_records


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("ids", help="text file containing ids of sequence to download from NCBI's nucleotide database (one name per line)")
    parser.add_argument("output", help="FASTA sequence corresponding to given ids")
    parser.add_argument("--accessions", action="store_true", help="input ids are accessions instead of nucleotide ids")
    args = parser.parse_args()

    with open(args.ids, "r") as fh:
        ids = [line.strip() for line in fh]

    # If the user supplied accessions, look up the nucleotide ids corresponding
    # to each accession.
    if args.accessions:
        nucleotide_ids = get_ids_by_accession(ids)
    else:
        nucleotide_ids = ids

    nucleotide_sequences = download_sequences(nucleotide_ids)
    SeqIO.write(nucleotide_sequences, args.output, "fasta")
