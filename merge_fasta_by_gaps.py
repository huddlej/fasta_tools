#!/bin/env python
"""
"""
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re

# Zero-based index of PSL array output corresponding to the size of the target
# sequence (i.e., the output sequence).
PSL_TSIZE_INDEX=10

def merge_fasta_by_gaps(input_filename, output_filename, sequence_name, gap_size):
    """
    Compress all consecutive instances of N bases from the given input FASTA file into a single N base.
    """
    gap = "N" * gap_size
    new_sequences = []
    chain = []
    output_position = 0
    for record in SeqIO.parse(input_filename, "fasta"):
        record_length = len(record)
        new_sequences.append(str(record.seq))

        # Create a simple PSL record for the input sequence and the new output
        # sequence for future liftOver use.
        chain.append([0, 0, 0, 0, 0, 0, 0, 0, "+", sequence_name, 0, output_position, output_position + record_length, record.id, record_length, 0, record_length, 1, record_length, 0, 0])

        output_position = output_position + record_length + gap_size

    # Save the new record with all sequences joined by the requested gaps.
    new_record = SeqRecord(Seq(gap.join(new_sequences)), id=sequence_name, description="")

    # Write out all records to the given output file.
    SeqIO.write([new_record], output_filename, "fasta")

    # Update all chain items with the final target sequence length and output the chain.
    new_record_length = len(new_record)
    for item in chain:
        item[PSL_TSIZE_INDEX] = new_record_length
        print "\t".join(map(str, item))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file")
    parser.add_argument("output_file")
    parser.add_argument("sequence_name", help="name to assign to output sequence")
    parser.add_argument("--gap_size", type=int, default=5000, help="size of gaps in bases to insert between each input FASTA")
    args = parser.parse_args()

    merge_fasta_by_gaps(args.input_file, args.output_file, args.sequence_name, args.gap_size)
