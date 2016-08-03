# Nik Krumm, 2012
# see documentation at : https://eichlerlab.gs.washington.edu/help/wiki/doku.php?id=software:fastq_to_bed.py

import sys
import os
import pysam
from subprocess import Popen, PIPE
import time
import csv
import operator


def gzip_read(folder_path, barcode):
	f = Popen('zcat ' + folder_path, shell=True, stdout=PIPE)
	read = pysam.AlignedRead()
	
	for line in f.stdout:
		read.qname = "%s:%s" % (barcode,line.strip("\n"))
		read.seq = f.stdout.next().strip("\n")
		_ = f.stdout.next()
		read.qual = f.stdout.next().strip("\n")
		yield read


with open(sys.argv[1],'r') as sample_file:
	r = csv.DictReader(sample_file,  delimiter='\t')
	input = {}
	for row in r:
		if row["sampleID"] not in input:
			input[row["sampleID"]] = {}
			input[row["sampleID"]]["out_bam"] = row["out_bam"]
			input[row["sampleID"]]["readgroups"] = []
		
		# add readgroup info
		if ("fastq_in_1" in row) and ("fastq_in_2" in row) and  (row["fastq_in_1"] != "") and (row["fastq_in_2"] != ""):
			input[row["sampleID"]]["readgroups"].append({"files": 2, "barcode": row["barcode"], "fastq_in_1": row["fastq_in_1"],"fastq_in_2": row["fastq_in_2"]})
			# do file checks
		else:
			input[row["sampleID"]]["readgroups"].append({"files": 1, "barcode": row["barcode"], "fastq_in_1": row["fastq_in_1"]})
			# do file checks

for sampleID in input:
	readgroups = [{'PU': rg["barcode"], 'ID': rg["barcode"], 'SM': sampleID,'CN': 'UW_GS_Eichler_Lab'} for rg in input[sampleID]["readgroups"]]
	
	header = {'RG': readgroups, 'SQ': [{'SN':'ref','LN':1}], 'HD':{'VN':1.0, 'SO': 'unsorted'},'PG':[{'ID':'fastq_to_bam'}], 'CO': ["This is an UNALIGNED bam file created from FASTQ files using fastq_to_bam.py","Read 1 and Read 2 are interleaved!", "Nik Krumm 2012 -- nkrumm@uw.edu"]}
	
	
	# open bamfile here
	outbam =  pysam.Samfile(input[sampleID]["out_bam"], "wb", header= header)
	
	t1 = time.time()
	cnt = 0
	for rg in input[sampleID]["readgroups"]:
		
		if rg["files"] ==2:
			source_1 = gzip_read(rg["fastq_in_1"], rg["barcode"])
			source_2 = gzip_read(rg["fastq_in_2"], rg["barcode"])
			while 1:
				try:
					read_1 = source_1.next()
					read_2 = source_2.next()
				except StopIteration:
					try: source_1.next()
					except:
						try: source_2.next()
						except: break		
						print "ERROR: source 2 still has unwritten reads! Read 1 and Read 2 did not have the same number of reads!"
					print "ERROR: source 1 still has unwritten reads! Read 1 and Read 2 did not have the same number of reads!"
				
				_ = outbam.write(read_1)
				_ = outbam.write(read_2)
				cnt +=2
		
		else:
			source_1 = gzip_read(rg["fastq_in_1"], rg["barcode"])
			source_2 = None
			while 1:
				try:
					read_1 = source_1.next()
				except StopIteration:
					break
				
				_ = outbam.write(read_1)
				cnt +=1
	
	t2 = time.time()
	print "Sample %s done in %f; total reads = %d" % (sampleID, t2-t1, cnt)
	
	outbam.close()
	
	del source_1
	del source_2
