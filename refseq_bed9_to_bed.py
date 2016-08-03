#!/bin/env python
"""
Parse the SQL data dump from UCSC's refGene.txt file into a BED file with one
line per exon.
"""
import argparse
import csv
import pprint
import sys

CHROM=2
STRAND=3
CODING_START=6
CODING_END=7
EXON_COUNT=8
EXON_STARTS=9
EXON_ENDS=10
SCORE=11
NAME=12


def main(filename, emit_introns):
    gene_count = 0

    with open(filename, "r") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            gene_count += 1
            coding_start = int(row[CODING_START])
            coding_end = int(row[CODING_END])
            exon_count = int(row[EXON_COUNT])
            exon_starts = map(int, row[EXON_STARTS].split(",")[:exon_count])
            exon_ends = map(int, row[EXON_ENDS].split(",")[:exon_count])

            if not emit_introns:
                exons = zip(exon_starts, exon_ends)
                for exon in exons:
                    # At the last coding exon and 3' UTR.
                    if exon[0] < coding_end and exon[1] > coding_end:
                        print "\t".join(map(str, [row[CHROM], exon[0], coding_end, "%s|%i" % (row[NAME], gene_count), row[SCORE], row[STRAND]]))
                        print "\t".join(map(str, [row[CHROM], coding_end, exon[1], "%s|%i" % (row[NAME], gene_count), row[SCORE], row[STRAND]]))
                    # At the 5' UTR and first coding exon
                    elif exon[0] < coding_start and exon[1] > coding_start:
                        print "\t".join(map(str, [row[CHROM], exon[0], coding_start, "%s|%i" % (row[NAME], gene_count), row[SCORE], row[STRAND]]))
                        print "\t".join(map(str, [row[CHROM], coding_start, exon[1], "%s|%i" % (row[NAME], gene_count), row[SCORE], row[STRAND]]))
                    else:
                        print "\t".join(map(str, [row[CHROM], exon[0], exon[1], "%s|%i" % (row[NAME], gene_count), row[SCORE], row[STRAND]]))
            else:
                # Emit introns which are the spaces between exon blocks (one
                # exon's end to the next exon's start).
                exons = zip(exon_ends[:-1], exon_starts[1:])
                for exon in exons:
                    print "\t".join(map(str, [row[CHROM], exon[0], exon[1], "%s|%i" % (row[NAME], gene_count), row[SCORE], row[STRAND]]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("refGene", help="refGene.txt from UCSC MySQL database dump")
    parser.add_argument("--introns", action="store_true", help="emit introns instead of exons")
    args = parser.parse_args()

    main(args.refGene, args.introns)
