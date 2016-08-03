"""
Given a multiple sequence alignment, calculate the identity for each pair of
sequences.
"""
import argparse
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("alignment", help="Multiple sequence alignment input in FASTA format")
    args = parser.parse_args()

    # Load the multiple sequence alignment and sort records by name in ascending
    # order.
    original_alignment = AlignIO.read(args.alignment, "fasta")
    alignment = MultipleSeqAlignment(sorted([record for record in original_alignment], key=lambda record: record.name))

    for sequence_j in alignment:
        for sequence_k in alignment:
            if sequence_j != sequence_k:
                total = 0
                matches = 0
                for i in xrange(alignment.get_alignment_length()):
                    total += 1

                    if sequence_j[i].upper() == sequence_k[i].upper():
                        matches += 1

                print "\t".join((sequence_j.name, sequence_k.name, str(matches / float(total))))
