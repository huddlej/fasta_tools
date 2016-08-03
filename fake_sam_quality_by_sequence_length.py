import argparse
import pysam


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_bam", help="BAM to fix")
    parser.add_argument("output_bam", help="filtered BAM")
    args = parser.parse_args()

    input_bam = pysam.AlignmentFile(args.input_bam, "rb")
    output_bam = pysam.AlignmentFile(args.output_bam, "wb", template=input_bam)

    for alignment in input_bam:
        alignment.qual = "~" * len(alignment.query_sequence)
        output_bam.write(alignment)

    input_bam.close()
    output_bam.close()
