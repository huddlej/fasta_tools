import argparse
import os
import pysam


def pair_reads_in_bam(input_samfile, output_samfile, unpaired_samfile):
    read_pairs = {}
    count = 0
    reads_paired = 0

    for read in input_samfile:
        if read.is_paired:
            if count % 10000 == 0:
                print "%i reads waiting for mate, %i paired" % (len(read_pairs.keys()), reads_paired)

            if len(read_pairs.keys()) > args.max_unpaired_reads:
                print "Dumping unpaired reads to disk"

                if unpaired_samfile is not None:
                    # Write unpaired reads to a separate BAM file for later sorting by samtools.
                    for read in read_pairs.values():
                        unpaired_samfile.write(read)

                # Reset the read pair index. All mates we haven't found by now
                # will get written to this same file later.
                del read_pairs
                read_pairs = {}

                # TODO: remove me!
                break

            if not read.qname in read_pairs:
                read_pairs[read.qname] = read
            else:
                output_samfile.write(read)
                output_samfile.write(read_pairs[read.qname])
                del read_pairs[read.qname]
                reads_paired += 1

            count += 1

    # Write remaining unpaired reads to a separate BAM file for later sorting by
    # samtools.
    for read in read_pairs.values():
        unpaired_samfile.write(read)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mapped_bam", help="Input BAM with reads in coordinate sorted order")
    parser.add_argument("paired_bam", help="Output BAM with reads in read sorted order")
    parser.add_argument("unpaired_bam", help="Temporary BAM to store unpaired reads")
    parser.add_argument("--max_unpaired_reads", type=int, default=500000, help="Maximum number of unpaired reads to keep in memory before writing to disk")
    args = parser.parse_args()

    samfile = pysam.AlignmentFile(args.mapped_bam, "rb")
    pairedreads = pysam.AlignmentFile(args.paired_bam, "wb", template=samfile)
    unpairedreads = pysam.AlignmentFile(args.unpaired_bam, "wb", template=samfile)

    # Pair reads for input BAM.
    pair_reads_in_bam(samfile, pairedreads, unpairedreads)

    # Sort unpaired reads by name.
    unpairedreads.close()
    sorted_unpaired_bam = args.unpaired_bam.replace(".bam", ".sorted.bam")
    pysam.sort("-n", args.unpaired_bam, sorted_unpaired_bam.rstrip(".bam"))
    sorted_unpaired_reads = pysam.AlignmentFile(sorted_unpaired_bam, "rb")

    # Pair reads in sorted temporary BAM without tracking unpaired reads again.
    unpairedreads = pysam.AlignmentFile(args.unpaired_bam, "wb", template=samfile)
    pair_reads_in_bam(sorted_unpaired_reads, pairedreads, unpairedreads)

    # Close files.
    sorted_unpaired_reads.close()
    unpairedreads.close()
    pairedreads.close()
    samfile.close()

    # Remove temporary file used for sorting.
    os.unlink(sorted_unpaired_bam)
