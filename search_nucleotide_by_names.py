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
Entrez.email = ""

# Set maximum records to return from NCBI's Entrez e-search to the maximum
# allowed by the service.
MAX_RECORDS_TO_RETURN = 50000


def search_clone_by_name(all_clone_names, return_accessions=False):
    """
    Search NCBI's clone database for the given clone names and return a list of
    tuples containing titles and accessions from the nucleotide database.

    >>> search_clone_by_name(["ABC13-48822500O14"])
    [('Homo sapiens FOSMID clone ABC13-48822500O14 from chromosome 4, complete sequence', 'AC216972')]
    >>> search_clone_by_name(["fakeclone"])
    []
    """
    # By default, expect no accession_clone_pairs.
    accession_clone_pairs = []

    # Process clones in batches to avoid overly long URL requests.
    total_clones = len(all_clone_names)
    step = 100
    clone_range = range(0, total_clones, step)

    for start in clone_range:
        end = min(start + step, total_clones)
        sys.stdout.write("Processing batch %i to %i\n" % (start, end))
        clone_names = all_clone_names[start:end]

        # Search NCBI clone database by clone name.
        try:
            handle = Entrez.esearch(db="clone", term=" OR ".join(clone_names), retmax=MAX_RECORDS_TO_RETURN)
            clone_record = Entrez.read(handle)
            sys.stderr.write("Found %s matching records\n" % clone_record["Count"])

            if int(clone_record["Count"]) > 0:
                sys.stderr.write("Finding links from clone database to nucleotide database\n")
                handle = Entrez.elink(dbfrom="clone", db="nuccore", id=",".join(clone_record["IdList"]))
                link_record = Entrez.read(handle)

                if len(link_record) > 0 and len(link_record[0]["LinkSetDb"]) > 0:
                    # Collect nucleotide database ids from eLink record.
                    nuccore_ids = [i["Id"] for i in link_record[0]["LinkSetDb"][0]["Link"]]

                    if return_accessions:
                        accession_clone_pairs.extend(nuccore_ids)
                    else:
                        # Get the summary of the nucleotide record(s) matching the id.
                        sys.stderr.write("Looking up nucleotide summaries\n")
                        handle = Entrez.esummary(db="nuccore", id=",".join(nuccore_ids), retmax=MAX_RECORDS_TO_RETURN)
                        nuccore_summary = Entrez.read(handle)
                        sys.stderr.write("Retrieved %s matching summaries\n" % len(nuccore_summary))

                        # Report the accession(s) for this term.
                        accession_clone_pairs.extend([(summary["Title"], summary["Caption"]) for summary in nuccore_summary])
        except urllib2.URLError:
            accession_clone_pairs.extend([("ncbi_lookup_failed", "")])

    return accession_clone_pairs


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
    parser.add_argument("clone_names", help="name of a clone to search for in NCBI")
    parser.add_argument("output", help="tab-delimited output of clone name/accession pairs OR sequence corresponding to clone names")
    parser.add_argument("--sequence", action="store_true", help="download sequence corresponding to clone names instead of returning accessions")
    args = parser.parse_args()

    with open(args.clone_names, "r") as fh:
        clone_names = [line.strip() for line in fh]

    # Query NCBI by clone names.
    if args.sequence:
        nucleotide_ids = search_clone_by_name(clone_names, return_accessions=True)
        nucleotide_sequences = download_sequences(nucleotide_ids)
        SeqIO.write(nucleotide_sequences, args.output, "fasta")
    else:
        accession_clone_pairs = search_clone_by_name(clone_names)

        # Find NCBI records that correspond to the requested clones. Entrez results
        # can return spurious results for clones that don't have sequence in the
        # database, so this step is a sanity check against returning false
        # accessions.
        with open(args.output, "w") as oh:
            writer = csv.writer(oh, delimiter="\t", lineterminator="\n")

            for clone in clone_names:
                accession_found = False
                for pair in accession_clone_pairs:
                    # If title contains the clone name, print the name/accession
                    # pair.
                    if clone in pair[0]:
                        accession_found = True
                        writer.writerow((clone, pair[1]))

                # Let the user know if the clone name didn't have an accession.
                if not accession_found:
                    writer.writerow((clone, "N/A"))
