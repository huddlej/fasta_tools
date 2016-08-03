#!/usr/bin/env python
"""
Fix order of columns in blast output to create "alignment" format output:

query start end length subject start end length other-fields...
"""
import csv
import sys

filename = sys.argv[1]
query_length = sys.argv[2]
subject_length_file = sys.argv[3]

fh = open(subject_length_file, "r")
reader = csv.reader(fh, delimiter="\t")
subject_length_by_name = dict([(row[0], row[1]) for row in reader])
fh.close()

fh = open(filename, "r")
csv_reader = csv.reader(fh, delimiter="\t")

QUERY=0
SUBJECT=1
QUERY_START=6
QUERY_END=7
SUBJECT_START=8
SUBJECT_END=9

for row in csv_reader:
    if row[QUERY].startswith("Query"):
        first_row = True
    else:
        first_row = False

    new_row = [
        row[QUERY],
        row[QUERY_START],
        row[QUERY_END],
        query_length if not first_row else "Query length",
        row[SUBJECT],
        row[SUBJECT_START],
        row[SUBJECT_END],
        subject_length_by_name.get(row[SUBJECT], "0") if not first_row else "Subject length"
    ]
    new_row = new_row + row[2:6] + row[10:]
    print "\t".join(new_row)
