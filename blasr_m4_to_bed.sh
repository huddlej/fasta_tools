#!/bin/bash

awk 'OFS="\t" { if ($5 != $9) { strand="-"; start=$12-$11; end=$12-$10 } else { strand="+"; start=$10; end=$11 } score=sprintf("%.1f\n", $4) * 10; print $2,start,end,$1,score,strand }' \
    | sed 's/\/0_[0-9]\+//'
