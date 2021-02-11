#!/bin/bash

BCFT=$1
FILE=$2
MDIR=$3/ld_easy_temp

## final filtering
echo "Starting final filtering."
$BCFT view -R $MDIR/regions.txt $FILE | bgzip -c > $MDIR/final.vcf.gz
echo "All done. Output file is:"
realpath $MDIR/final.vcf.gz
