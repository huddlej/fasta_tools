import argparse
from collections import defaultdict
import pysam


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("barcodes", help="FASTQ of barcodes to count")
    args = parser.parse_args()

    fastq = pysam.FastqFile(args.barcodes)
    counts_by_barcode = defaultdict(int)

    for record in fastq:
        counts_by_barcode[record.sequence] += 1

    for key, value in counts_by_barcode.iteritems():
        print "%s\t%s" % (key, value)
