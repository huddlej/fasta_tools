#!/bin/env python
"""
"""
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import csv
import pysam
import re

# Zero-based index of PSL array output corresponding to the size of the target
# sequence (i.e., the output sequence).
PSL_TSIZE_INDEX=10


def write_records(sequence_name, new_sequences, output, chain):
    # Save the new record with all sequences joined by the requested gaps and
    # write out all records to the given output file.
    new_record = SeqRecord(Seq("".join(new_sequences).rstrip("N")), id=sequence_name, description="")
    SeqIO.write((new_record,), output, "fasta")

    # Update all chain items with the final target sequence length and output the chain.
    new_record_length = len(new_record)
    for item in chain:
        item[PSL_TSIZE_INDEX] = new_record_length
        print "\t".join(map(str, item))


def merge_fasta_by_gaps(assemblies, tiling_path, output_filename, chromosome_suffix, default_gap_size):
    """
    Merge all FASTA records from the same organismal chromosome into a single
    sequence.
    """
    current_chromosome = None
    assemblies = pysam.FastaFile(assemblies)
    tiling_path_fh = open(tiling_path, "r")
    tiling_path_reader = csv.reader(tiling_path_fh, delimiter="\t")
    output = open(output_filename, "w")

    for path in tiling_path_reader:
        chromosome, start, end, contig, contig_start, contig_end = path
        start, end, contig_start, contig_end = map(int, (start, end, contig_start, contig_end))
        contig_length = assemblies.get_reference_length(contig)
        path_length = contig_end - contig_start

        if current_chromosome != chromosome:
            if current_chromosome is not None:
                write_records(new_chromosome, new_sequences, output, chain)

            new_chromosome = "%s%s" % (chromosome, chromosome_suffix)
            new_sequences = []
            chain = []
            output_position = 0
            current_chromosome = chromosome
            previous_end = None

        # If the end of the last interval is equal to the start of this
        # interval, set the gap size to zero. Otherwise use the default gap
        # size.
        if previous_end == start:
            gap_size = 0
        else:
            gap_size = default_gap_size

        if path_length > 0:
            new_sequences.append(assemblies.fetch(contig, contig_start, contig_end) + gap_size * "N")
            chain.append([0, 0, 0, 0, 0, 0, 0, 0, "+", new_chromosome, None, output_position, output_position + path_length, contig, contig_length, contig_start, contig_end, 1, path_length, 0, 0])
            output_position = output_position + path_length + gap_size

        # Set the previous end for the next iteration.
        previous_end = end

    # Write out the final batch of sequences.
    if len(new_sequences) > 0:
        write_records(new_chromosome, new_sequences, output, chain)

    tiling_path_fh.close()
    output.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("assemblies", help="FASTA file of assemblies with a tiling path")
    parser.add_argument("tiling_path", help="tiling path of assemblies in BED 3+3 format with reference coordinates in the first three columns and corresponding assembly coordinates in the last three columns")
    parser.add_argument("output_file", help="FASTA file of assemblies merged by gaps across the given tiling path")
    parser.add_argument("chromosome_suffix", help="suffix to add to each chromosome name in the output sequence to distinguish it from the standard reference")
    parser.add_argument("--gap_size", type=int, default=1000, help="size of gaps in bases to insert between each input FASTA")
    args = parser.parse_args()

    merge_fasta_by_gaps(args.assemblies, args.tiling_path, args.output_file, args.chromosome_suffix, args.gap_size)
