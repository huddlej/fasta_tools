#!/usr/bin/env python

import sys

number_of_ns = int(sys.argv[1])
width = int(sys.argv[2])
remainder = number_of_ns % width

for n in xrange(number_of_ns / width):
    print "N"*width

if remainder > 0:
    print "N"*remainder
