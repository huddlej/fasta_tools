"""
Convert a raw .qual database into a FASTA file using FASTQ encoding of quality
scores.

>CH17-224P20.FORWARD.1
10  6  5  5  6  7  7  7  4  5 10 11  6  5  6  8  8  8  8  8
 7  9  7  8  9 21 10 16 24 11  7 18  8  7  8  9 24 34 33 38
11 10 21  7 13  8  9 11 19 13 23 19 24 59 41 52 13 20 15 29
"""
import argparse
import cPickle
import re

FASTQ_QUALITIES = """!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI"""
MAX_FASTQ_QUALITY = 40

def fastq2qual(fastq_string):
    """
    Converts the given FASTQ quality string to a tuple of quality integers.

    >>> fastq2qual('!"#II!')
    (0, 1, 2, 40, 40, 0)
    """
    return tuple([FASTQ_QUALITIES.index(character) for character in fastq_string])


def qual2fastq(qualities):
    """
    Converts the given tuple of quality integers into a FASTQ quality string.

    >>> qual2fastq((0, 1, 2, 40, 40, 0))
    '!"#II!'
    """
    return "".join([FASTQ_QUALITIES[min(quality, MAX_FASTQ_QUALITY)] for quality in qualities])


def quality_database_to_fastq_bed(quality_database):
    sequence_name = None

    with open(quality_database, "r") as quality_fh:
        for line in quality_fh:
            if line.startswith(">"):
                if sequence_name is not None:
                    print "\t".join(map(str, (sequence_name, 1, len(sequence), sequence)))

                sequence_name = line.rstrip()[1:]
                sequence = ""
            else:
                sequence += qual2fastq([int(i) for i in re.split("\s+", line.strip())])

        print "\t".join(map(str, (sequence_name, 1, len(sequence), sequence)))


def qual2pickle(quality_database, pickle_file):
    sequence_name = None
    qualities_by_clone_end = {}

    with open(quality_database, "r") as quality_fh:
        for line in quality_fh:
            if line.startswith(">"):
                if sequence_name is not None:
                    qualities_by_clone_end[sequence_name] = tuple(qualities_by_clone_end[sequence_name])

                sequence_name = line.rstrip()[1:]
                qualities_by_clone_end[sequence_name] = []
            else:
                qualities_by_clone_end[sequence_name].extend([int(i) for i in re.split("\s+", line.strip())])

    qualities_by_clone_end[sequence_name] = tuple(qualities_by_clone_end[sequence_name])

    with open(pickle_file, "wb") as pickle_fh:
        cPickle.dump(qualities_by_clone_end, pickle_fh)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("command")
    parser.add_argument("input", nargs="+")
    args = parser.parse_args()

    if args.command == "qual2fastq":
        quality_database_to_fastq_bed(args.input[0])
    elif args.command == "qual2pickle":
        qual2pickle(*args.input)
