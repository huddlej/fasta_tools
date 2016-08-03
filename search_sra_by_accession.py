#!/bin/env python
"""
Given a clone name, report any matching accessions in NCBI's clone database.
"""
import argparse
from Bio import Entrez
import urllib2
from xml.etree import ElementTree
Entrez.email = ""

MAX_RECORDS=1000


def search_sra_by_accessions(accessions):
    """
    Search NCBI's SRA database for the given experiment accessions and return a list of
    tuples with experiment and run accessions.

    >>> search_sra_by_accessions(["SRX069099"])
    [('SRX069099', 'SRR231127')]
    >>> search_sra_by_accessions(["fakeaccession"])
    []
    """
    # By default, expect no accessions.
    run_accessions = []

    # Search SRA database by experiment accessions.
    search_accessions = " OR ".join(accessions)

    try:
        handle = Entrez.esearch(db="sra", term=search_accessions, retmax=MAX_RECORDS)
        search_records = Entrez.read(handle)

        if search_records["Count"] > 0:
            experiment_ids = ",".join(search_records["IdList"])
            handle = Entrez.efetch(db="sra", id=experiment_ids, retmax=MAX_RECORDS)
            experiment_xml_string = "\n".join(handle.readlines())
            experiment_xml = ElementTree.fromstring(experiment_xml_string)

            experiment_packages = experiment_xml.findall("EXPERIMENT_PACKAGE")
            for experiment_package in experiment_packages:
                experiment = experiment_package.find("EXPERIMENT")
                experiment_accession = experiment.attrib["accession"]
                run_accessions.extend([(experiment_accession, run.attrib["accession"]) for run in experiment_package.findall("RUN_SET/RUN")])
    except urllib2.URLError:
        run_accessions = ["ncbi_lookup_failed"]

    return run_accessions


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("accessions_file", help="the name of a file containing one or more SRA experiment accessions (e.g., SRX069099) per line")
    args = parser.parse_args()

    # Load experiment accessions from the given file.
    with open(args.accessions_file, "r") as accessions_fh:
        experiment_accessions = [line.strip() for line in accessions_fh]

    accessions = search_sra_by_accessions(experiment_accessions)
    print "\n".join(["\t".join(accession_pair) for accession_pair in accessions])
