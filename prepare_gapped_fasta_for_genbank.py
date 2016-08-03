#!/bin/env python
"""
Prepare a FASTA with gaps represented as Ns for submission to GenBank.
"""
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re


def prepare_gapped_fasta_for_genbank(input_filename, output_filename):
    # Load the original FASTA sequence.
    fasta = SeqIO.parse(input_filename, "fasta")
    original_record = fasta.next()
    original_sequence = str(original_record.seq)

    # Split the FASTA sequence on one or more Ns to get individual segments.
    segments = re.split("N+", original_sequence)

    # Prepare a new named record for the first segment.
    new_segments = [SeqRecord(Seq(segments[0]), id=original_record.id, description="")]

    # Create gap records for the remaining segments.
    segment_count = 2
    for segment in segments[1:]:
        new_segments.append(SeqRecord(Seq(segment), id="segment%i" % segment_count, description=""))
        segment_count += 1

    # Write out a new FASTA file for all segments.
    SeqIO.write(new_segments, output_filename, "fasta")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file")
    parser.add_argument("output_file")
    args = parser.parse_args()

    prepare_gapped_fasta_for_genbank(args.input_file, args.output_file)
