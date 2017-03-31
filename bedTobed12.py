#!/bin/env python
import argparse
import csv
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# create console handler and set level to debug
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)

# create formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# add formatter to ch
ch.setFormatter(formatter)

# add ch to logger
logger.addHandler(ch)

CHROMOSOME=0
START=1
END=2
FEATURE=3
STRAND=5


def bed_to_bed12(filename, ignore_strand):
    """
    Convert a BED file to a BED12 file for each feature in the BED file across
    an entire chromosome.

    Assumes data are sorted by chromosome, feature, and start position.
    """
    chromosome = None
    feature = None

    # Group coordinates by chromosome and then by feature.
    blocks_by_chromosome = {}
    with open(filename, "r") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if not ignore_strand and len(row) > STRAND - 1:
                strand = row[STRAND]
            else:
                strand = "+"

            feature_strand=";".join((row[FEATURE], strand))
            blocks_by_chromosome.setdefault(row[CHROMOSOME], {}).setdefault(feature_strand, []).append(list(map(int, row[START:END+1])))

    # Print a line for each set of blocks per chromosome/feature combination.
    score = 0
    color = "255,0,0"
    for chromosome, features in blocks_by_chromosome.items():
        for feature_strand, blocks in features.items():
            feature, strand = feature_strand.split(";")
            flat_blocks = [item for sublist in blocks for item in sublist]
            start = min(flat_blocks)
            end = max(flat_blocks)
            block_sizes = [block[1] - block[0] for block in blocks]
            block_starts = [block[0] - start for block in blocks]

            # chr8    2697394 2698108 172343_ABC9_3_5_000043351700_N11        0       +       2697394 2698108 255,0,0 2       714,683 0,29
            if (start + block_starts[-1] + block_sizes[-1]) != end:
                block_sizes = [end - start]
                block_starts = [0]

            print("\t".join(map(str, [chromosome, start, end, feature, score, strand, start, end, color, len(block_sizes), ",".join(map(str, block_sizes)), ",".join(map(str, block_starts))])))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_bed")
    parser.add_argument("--ignore_strand", action="store_true")
    args = parser.parse_args()

    bed_to_bed12(args.input_bed, ignore_strand=args.ignore_strand)
