#!/bin/env python
"""
Convert genomicSuperDup files from WGAC into standard GFF3 format for use by
NCBI in their browser.
"""
import argparse
import csv

# Define constants.
FIELDS = [
    "chrom",
    "chromStart",
    "chromEnd",
    "name",
    "score",
    "strand",
    "otherChrom",
    "otherStart",
    "otherEnd",
    "otherSize",
    "uid",
    "posBasesHit",
    "testResult",
    "verdict",
    "chits",
    "ccov",
    "alignfile",
    "alignL",
    "indelN",
    "indelS",
    "alignB",
    "matchB",
    "mismatchB",
    "transitionsB",
    "transversionsB",
    "fracMatch",
    "fracMatchIndel",
    "jcK",
    "k2K"
]

def convert_wgac_row_to_gff3(row, row_id):
    """
    Given a row of WGAC data, convert it to a GFF3 entry.
    """
    row_dict = dict(zip(FIELDS, row))

    # Convert zero-based coordinates to one-based coordinates.
    row_dict["chromStart"] = str(int(row_dict["chromStart"]) + 1)
    row_dict["otherStart"] = str(int(row_dict["otherStart"]) + 1)

    gff3_row = [
        row_dict["chrom"],
        "wgac",
        "match",
        row_dict["chromStart"],
        row_dict["chromEnd"],
        "0",
        row_dict["strand"].replace("_", "-"),
        "."
    ]

    # Build GFF3-specific attributes.
    row_dict["ID"] = str(row_id)
    row_dict["Name"] = row_dict["name"]
    row_dict["Target"] = " ".join((
        row_dict["otherChrom"],
        row_dict["otherStart"],
        row_dict["otherEnd"]
    ))

    # Construct the attributes field from all GFF3 attributes plus all other
    # WGAC attributes.
    attributes = ";".join(
        ["=".join((field, row_dict[field]))
         for field in ["ID", "Name", "Target"] + FIELDS]
    )
    gff3_row.append(attributes)

    return gff3_row


def convert_wgac_to_gff3(wgac_filename):
    with open(wgac_filename, "r") as fh:
        reader = csv.reader(fh, delimiter="\t")

        # Emit the GFF3 header.
        print "##gff-version 3"

        row_id = 0
        for row in reader:
            print "\t".join(convert_wgac_row_to_gff3(row, row_id))
            row_id += 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("wgac", help="genomicSuperDups file from WGAC pipeline")
    args = parser.parse_args()

    convert_wgac_to_gff3(args.wgac)
