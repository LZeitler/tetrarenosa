#!/bin/bash

BCFT=$1
FILE=$2
SCAF=$3
MDIR=$4/ld_easy_temp

mkdir -p $MDIR

## filtering and make matrix
echo "Starting filtering."
$BCFT view -r $SCAF -m 2 -M 2 -q 0.05 -Q 0.95 -U -i 'F_MISSING<0.1' $FILE | bgzip -c > $MDIR/chrom_$SCAF.vcf.gz
if [ ! -f ${MDIR}/names.txt ]; then
    zgrep "#CHROM" $MDIR/chrom_$SCAF.vcf.gz | awk '{print substr($0, index($0,$10))}' > $MDIR/names.txt
    echo "Individual names written."
fi
$BCFT annotate -x INFO,^FORMAT/GT $MDIR/chrom_$SCAF.vcf.gz | grep -v '#' | bgzip -c > $MDIR/genos_$SCAF.vcf.gz
echo "Writing genotype matrix and positions."
zcat $MDIR/genos_$SCAF.vcf.gz | awk '{print $2}' | less -S > $MDIR/positions_$SCAF.txt
zcat $MDIR/genos_$SCAF.vcf.gz | awk '{print substr($0, index($0,$10))}' > $MDIR/matrix_$SCAF.txt

echo "Starting conversion."
sed -i -e 's;0/0/0/0;0;g'  $MDIR/matrix_$SCAF.txt
sed -i -e 's;0/0/0/1;1;g'  $MDIR/matrix_$SCAF.txt
sed -i -e 's;0/0/1/1;2;g'  $MDIR/matrix_$SCAF.txt
sed -i -e 's;0/1/1/1;3;g'  $MDIR/matrix_$SCAF.txt
sed -i -e 's;1/1/1/1;4;g'  $MDIR/matrix_$SCAF.txt
sed -i -e 's;./././.;NA;g' $MDIR/matrix_$SCAF.txt

echo "$SCAF done."
