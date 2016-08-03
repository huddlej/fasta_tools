#!/bin/env python
"""
Convert a MAST XML file (from the MEME suite) to a BED file that includes the
positions of each motif in each sequence.
"""
import argparse
import xml.etree.ElementTree as ET


def mast_to_bed(mast_xml):
    """
    Given the filename of a MAST XML file, print the corresponding BED output to
    standard out.
    """
    tree = ET.parse(mast_xml)
    root = tree.getroot()
    sequences = root.find("sequences").findall("sequence")

    for sequence in sequences:
        name = sequence.attrib["name"]

        for seg in sequence.findall("seg"):
           for hit in seg.findall("hit"):
               start = hit.attrib["pos"]
               end = str(int(start) + len(hit.attrib["match"]))
               motif = hit.attrib["motif"]

               print "\t".join((name, start, end, motif))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mast_xml")
    args = parser.parse_args()

    mast_to_bed(args.mast_xml)
