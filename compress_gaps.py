#!/bin/env python
"""
"""
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re


def compress_gaps(input_filename, output_filename, trim_gaps):
    """
    Compress all consecutive instances of N bases from the given input FASTA file into a single N base.
    """
    if trim_gaps:
        gap_replacement = ""
    else:
        gap_replacement = "N"

    with open(output_filename, "w") as oh:
        for record in SeqIO.parse(input_filename, "fasta"):
            # Remove all lowercase letters in the sequence for this record.
            new_sequence = re.sub("[nN]+", gap_replacement, str(record.seq))

            # If any sequence remains, create a new record for it.
            if len(new_sequence) > 0:
                # Write out all records to the given output file.
                SeqIO.write((SeqRecord(Seq(new_sequence), id=record.id, description=""),), oh, "fasta")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file")
    parser.add_argument("output_file")
    parser.add_argument("--trim", action="store_true", help="trim gaps completely instead of compressing to single bases")
    args = parser.parse_args()

    compress_gaps(args.input_file, args.output_file, args.trim)
