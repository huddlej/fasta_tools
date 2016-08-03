#!/bin/env python
from xml.etree import cElementTree as et
import sys
import urllib


def get_xml_by_clone_ids(clone_ids_file, output_dir):
    query = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nucleotide&id=%s"
    with open(clone_ids_file, "r") as clone_ids_fh:
        for clone_id in clone_ids_fh:
            clone_id = clone_id.strip()
            resource = urllib.urlopen(query % clone_id)

            with open("%s/%s.xml" % (output_dir, clone_id), "w") as clone_xml_fh:
                for line in resource:
                    clone_xml_fh.write(line)


def xml_to_tab(xml_list_filename):
    columns_to_print = ("Gi", "Caption", "Title", "Extra")
    with open(xml_list_filename, "r") as xml_list_fh:
        for xml_filename in xml_list_fh:
            with open(xml_filename.strip(), "r") as xml_fh:
                tree = et.parse(xml_fh)
                columns = dict([(item.get("Name"), item.text)
                                for item in tree.getroot().findall("DocSum/Item")
                                if item.get("Name") in columns_to_print])
                print "\t".join([columns.get(column, "") for column in columns_to_print])


if __name__ == "__main__":
    if len(sys.argv) == 2:
        # Convert XML to tab.
        xml_to_tab(sys.argv[1])
    elif len(sys.argv) == 3:
        # Get XML.
        get_xml_by_clone_ids(sys.argv[1], sys.argv[2])
    else:
        sys.stderr.write("Usage: %prog clone_ids.txt xml_output\n\n")
        sys.exit(1)


