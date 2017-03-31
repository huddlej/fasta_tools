#!/bin/bash
head -n 1 | awk 'BEGIN { OFS="\t"; FS="\t" } { for (i = 1; i <= NF; i++) { print i,$i } }'