#!/bin/bash

base_name=`basename $1`
/net/eichler/vol2/local/inhousebin/fasta_length -o /tmp/${base_name}.tmp.log $1 &> /dev/null
sed '1d;$d' /tmp/${base_name}.tmp.log | cut -f 3-4
rm -f /tmp/${base_name}.tmp.log
