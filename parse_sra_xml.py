#!/bin/env python

import csv

library = "CT3576"

oh = open("%s.manifest.tab" % library, "w")
writer = csv.writer(oh, delimiter="\t")

fh = open("%s.xml" % library, "r")
t = ElementTree.parse(fh)
r = t.getroot()

for run in r.findall("RUN_SET/RUN"):
    accession = run.attrib.get("accession")
    alias = run.attrib.get("alias")
    flowcell, lane = alias.split("_")[3].split("-")
    for run_number in range(1, 3):
        fastq = "/net/eichler/vol18/genomes/high_coverage_trios/NA18507_trio/NA18507/sequence/%s/fastq/%s_%s.fastq.gz" % (library, accession, run_number)
        row = (library, flowcell, lane, str(run_number), fastq)
        writer.writerow(row)

oh.close()
fh.close()
