#!/bin/env python
from xml.etree import cElementTree as et
import sys
import urllib


def get_accessions_by_clone_names(clone_names_file):
    query = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucleotide&term=%s"
    fh = open(clone_names_file, "r")

    for clone_name in fh:
        clone_name = clone_name.strip()
        resource = urllib.urlopen(query % clone_name)
        tree = et.parse(resource)
        root = tree.getroot()
        count = int(root.findall("Count")[0].text)
        id_list = root.findall("IdList")
        errors = root.findall("ErrorList")

        # Handle three different cases:
        # 1. Query failed (e.g., because clone doesn't exist) and an error was returned.
        # 2. Query returned more than one result and the user needs to choose one.
        # 3. Query returned exactly one result.
        if len(errors) > 0:
            sys.stderr.write("%s\tquery_error\n" % clone_name)
        elif len(id_list) > 0 and len(id_list[0].getchildren()) > 1:
            sys.stderr.write("%s\ttoo_many_results (%i)\n" % (clone_name, len(id_list[0].getchildren())))
        else:
            print "\t".join((clone_name, id_list[0].getchildren()[0].text))

        resource.close()

    fh.close


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.stderr.write("Usage: ./get_accessions_by_clone_names clonenames.txt > matches.tab 2> missing_clones.tab\n\n")
        sys.stderr.write(" - Matching accession ids are written to stdout.\n")
        sys.stderr.write(" - Clones without matches are written to stderr.\n")
        sys.exit(1)

    get_accessions_by_clone_names(sys.argv[1])
