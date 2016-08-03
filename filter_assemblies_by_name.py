import argparse
import pysam
import sys


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("assemblies")
    parser.add_argument("names")
    args = parser.parse_args()

    with open(args.names, "r") as fh:
        assembly_names = set([line.strip() for line in fh])

    sys.stderr.write("Found %i assembly names\n" % len(assembly_names))

    bam = pysam.AlignmentFile(args.assemblies, "rb")
    assemblies_kept = 0
    for alignment in bam:
        if alignment.query_name in assembly_names:
            assemblies_kept += 1
            print(">%s" % alignment.query_name)
            print(alignment.query_sequence)

    sys.stderr.write("Kept %i assembly sequences\n" % assemblies_kept)
