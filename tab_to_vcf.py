#!/bin/env python
"""
Convert a tab-delimited file into a valid VCF.
"""
import argparse
import csv
import pysam
import re


def convert_tab_to_vcf(data_filename, reference_filename, vcf_header_filename):
    # Load the reference.
    reference = pysam.Fastafile(reference_filename)

    # Print VCF header.
    with open(vcf_header_filename, "r") as fh:
        for line in fh:
            print line.rstrip()

    # Build initial column header.
    columns = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]

    with open(data_filename, "r") as fh:
        reader = csv.reader(fh, delimiter="\t")
        count = 0

        # Load species and header information.
        species_row = reader.next()
        header = reader.next()

        # Compute sample names with species names.
        sample_start_index = None
        species_sample_names = []
        for i in xrange(len(species_row)):
            if species_row[i] != "":
                if sample_start_index is None:
                    sample_start_index = i

                species_name = re.sub("\s+", "_", species_row[i].rstrip())
                species_sample_names.append("-".join((species_name, header[i])))

        # Add computed sample names to VCF header columns and emit the header.
        columns.extend(species_sample_names)
        print "\t".join(columns)

        # Setup unique id for autoincrementing.
        unique_id = 0

        # Set constants.
        quality = "."
        filter_pass = "PASS"
        format = "GT"

        for row in reader:
            # Setup fields from input.
            chromosome = row[0]
            start = row[1]
            end = row[2]
            record_type = row[3]
            record_id = "%s_%s" % (record_type, unique_id)
            reference_allele = reference.fetch(chromosome, int(start) - 1, int(start)).upper()

            if "L1" in record_type:
                alternate_allele = "<INS:ME:LINE1>"
            elif "Alu" in record_type:
                alternate_allele = "<INS:ME:ALU>"

            max_support = row[4]
            strand = row[5]
            total_support = row[6]
            length = row[7]

            # Build INFO field.
            info = {
                "SVLEN": length,
                "END": end,
                "SVTYPE": record_type,
                "MEINFO": ",".join((record_type, "1", str(length), strand)),
                "MAXSUPPORT": max_support,
                "TOTALSUPPORT": total_support
            }
            info = ";".join(["=".join(item) for item in info.iteritems()])

            # Build genotypes assigning a heterozygous genotype when sample has
            # a value of 1 and reference allele when sample is 0.
            genotypes = [int(sample) and "1/0" or "0/0" for sample in row[sample_start_index:]]

            # Emit VCF record.
            print "\t".join([chromosome, start, record_id, reference_allele, alternate_allele, quality, filter_pass, info, format] + genotypes)

            unique_id += 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("data")
    parser.add_argument("reference")
    parser.add_argument("vcf_header")
    args = parser.parse_args()

    convert_tab_to_vcf(args.data, args.reference, args.vcf_header)
