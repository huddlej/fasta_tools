import argparse
from Bio import SeqIO


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta")
    args = parser.parse_args()

    for name, record in SeqIO.index(args.fasta, "fasta").iteritems():
        print "\t".join((name, str(len(record))))

