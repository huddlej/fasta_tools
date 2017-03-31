#!/bin/bash

awk 'OFS="\t" { if ($3 != $4) { strand="-"; start=$9-$8; end=$9-$7 } else { strand="+"; start=$7; end=$8 } score=sprintf("%.1f\n", $6) * 10; print $2,start,end,$1,score,strand }' \
    | sed 's/\/0_[0-9]\+//'
