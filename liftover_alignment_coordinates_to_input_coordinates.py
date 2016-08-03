import argparse
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import csv
import pybedtools


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_alignment", help="Multiple sequence alignment in FASTA format")
    parser.add_argument("coordinates", help="BED file of coordinates in the space of the given alignment; the chromosome field is required but not used.")
    args = parser.parse_args()

    alignments = AlignIO.read(args.input_alignment, "fasta")
    coordinates = pybedtools.BedTool(args.coordinates)

    for alignment in alignments:
        # Find total gap bases preceding each position in the alignment. The
        # first position by definition has zero gaps in front of it.
        gaps = [0]
        for i in xrange(alignments.get_alignment_length()):
            if alignment[i] == "-":
                gaps.append(gaps[-1] + 1)
            else:
                gaps.append(gaps[-1])

        # Emit the given coordinates translated into the coordinate space of the
        # original input sequence (i.e., coordinates minus total gap bases
        # preceding each position).
        for coordinate in coordinates:
            print "\t".join((alignment.name, str(coordinate.start - gaps[coordinate.start]), str(coordinate.end - gaps[coordinate.end])))
