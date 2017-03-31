#!/bin/bash

# chr10_131587030_131607130_ctg7180000000002_rc/0_43248   43248   0       43248   +       chr10:131557030-131637130       80100   13692   61628   -       -190191 43034
# chr10_131587030_131607130_ctg7180000000002/0_43248      43248   0       43248   +       chr10:131557030-131637130       80100   13692   61628   +       -190191 43034

awk 'OFS="\t" {
    matches=0;
    mismatches=0;
    alignment=$18
    last_was_mismatch=0;

    for (i=1; i <= length(alignment); i++) {
        if (substr(alignment, i, 1) == "|") {
            matches += 1;
            last_was_mismatch=0;
        }
        else {
            if (last_was_mismatch == 0) {
                mismatches += 1
            }

            last_was_mismatch=1
        }
    }

    score=sprintf("%.3f", matches / (matches + mismatches)) * 1000;

    if ($5 != $10) {
        strand = "-";
        start = $7 - $9;
        end = $7 - $8
    }
    else {
        strand = "+";
        start = $8;
        end = $9
    }

    num_of_pieces=split($1, query_pieces, "/");
    query_name = query_pieces[1];
    subject_name = $6;
    print subject_name,start,end,query_name,score,strand
}'
