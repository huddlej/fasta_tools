import argparse
import pysam


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("bam")
    parser.add_argument("fasta")
    args = parser.parse_args()

    bam = pysam.AlignmentFile(args.bam, "rb")
    fasta = pysam.FastaFile(args.fasta)
    fixed_bam = pysam.AlignmentFile("-", "wb", template=bam)

    for alignment in bam:
        alignment.query_sequence = fasta.fetch(alignment.query_name)
        alignment.template_length = len(alignment.query_sequence)
        fixed_bam.write(alignment)
