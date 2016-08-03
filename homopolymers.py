#!/bin/env python
"""
Print coordinates of continuous matching bases.
"""
import fileinput

fasta = fileinput.input()
name=fasta.readline().strip().replace(">", "")
previous_base=""
base_count=1
start_position=0
current_position=0

for line in fasta:
    for base in line.strip():
        if previous_base == base or previous_base == "":
            base_count += 1
            previous_base = base
        else:
            print "\t".join(map(str, (name,start_position,current_position,previous_base)))
            start_position = current_position
            previous_base = base

        current_position += 1
