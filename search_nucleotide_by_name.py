#!/bin/env python
"""
Given a clone name, report any matching accessions in NCBI's nucleotide
database.
"""
import argparse
from Bio import Entrez
import urllib2
Entrez.email = ""


def search_clone_by_name(clone_name):
    """
    Search NCBI's clone database for the given clone name and return a list of
    zero or more accessions from the nucleotide database.

    >>> search_clone_by_name("1205403_ABC13_11_000048822500_O14")
    ['AC216972']
    >>> search_clone_by_name("fakeclone")
    []
    """
    # By default, expect no accessions.
    accessions = []

    # Search NCBI clone database by clone name.
    try:
        handle = Entrez.esearch(db="nuccore", term=clone_name)

        clone_record = Entrez.read(handle)

        if clone_record["Count"] > 0:
            # Get the summary of the nucleotide record(s) matching the clone id.
            handle = Entrez.esummary(db="nuccore", id=",".join(clone_record["IdList"]))
            nuccore_summary = Entrez.read(handle)

            # Report the accession(s) for this term.
            accessions = [summary["Caption"] for summary in nuccore_summary]
    except urllib2.URLError:
        accessions = ["ncbi_lookup_failed"]

    return accessions


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("clone_name", help="name of a clone to search for in NCBI")
    args = parser.parse_args()

    accessions = search_clone_by_name(args.clone_name)
    print ",".join(accessions)
