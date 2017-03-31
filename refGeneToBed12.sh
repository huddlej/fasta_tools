#!/bin/bash

awk 'OFS="\t" { exon_count=split($10, starts, ","); exon_count=split($11, ends, ","); sizes=ends[1] - starts[1]; new_starts=starts[1]-$5; for (i = 2; i <= $9; i++) { sizes = sizes","(ends[i]-starts[i]); new_starts=new_starts","starts[i]-$5 } print $3,$5,$6,$13,$12,$4,$7,$8,"0,0,0",$9,sizes,new_starts }'
