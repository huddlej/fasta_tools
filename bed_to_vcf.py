"""
chr1    66277   66377   deletion        100     TatattatataatatataatataaatataatataaattatataaatataatatatattttattatataatataatatatattatataaatataatataTA    (TATAA)n:FULL   no_tsd  chr1-20000-80000|ctg7180000000002|quiver/0_44021     26527   26627   0.97
chr1    90068   90127   insertion       59      AGACAGTCCCTCAGTCCCTCTGTCTCTGCCAATCAGTTAACCTGCTGCTTCCTGGAGGA     NONE    AGACAGTCCCTCAGTCCCT     chr1-40000-100000|ctg7180000000001|quiver/0_62902       46220   462790.00
"""
import argparse
import csv
import pysam


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("bed")
    parser.add_argument("reference")
    parser.add_argument("vcf")
    args = parser.parse_args()

    SMRTSV_TYPE_TO_SV_TYPE = {
        "insertion": "INS",
        "deletion": "DEL"
    }

    reference = pysam.FastaFile(args.reference)

    with open(args.vcf, "w") as vcf_fh:
        vcf_fh.write("""##fileformat=VCFv4.2\n""")
        vcf_fh.write("""##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n""")
        vcf = csv.writer(vcf_fh, delimiter="\t", lineterminator="\n")
        vcf.writerow(("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"))

        with open(args.bed, "r") as bed_fh:
            reader = csv.reader(bed_fh, delimiter="\t")
            for row in reader:
                CHROM = row[0]
                POS = row[1]
                ID = "1"
                QUAL = "30"
                FILTER = "PASS"
                INFO = [
                    ("END", row[2]),
                    ("SVTYPE", SMRTSV_TYPE_TO_SV_TYPE[row[3]]),
                    ("SVLEN", row[4]),
                    ("REPEATS", row[6]),
                    ("TSD", row[7]),
                    ("CONTIG", row[8]),
                    ("CONTIG_START", row[9]),
                    ("CONTIG_END", row[10])
                ]

                preceding_base = reference.fetch(reference=CHROM, start=int(POS), end=int(POS) + 1)

                if row[3] == "insertion":
                    REF = preceding_base
                    ALT = row[5]
                elif row[3] == "deletion":
                    REF = row[5]
                    ALT = preceding_base

                INFO = ";".join(["=".join(item) for item in INFO])
                vcf.writerow((CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO))

    reference.close()
