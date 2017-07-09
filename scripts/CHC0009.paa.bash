#! /bin/bash

PATH=/mnt/wigclust5/data/safe/kendall/bowtie-0.12.7/:/mnt/wigclust5/data/safe/kendall/samtools-0.1.16/:/mnt/wigclust5/data/safe/kendall/Python-2.7.1/:$PATH

BYTES_MEM=$(ls -l /mnt/wigclust5/data/safe/kendall/riggs01/barcode.split/FC64VNGAAXX/s_5_sequence.WB1.txt.gz | cut -d " " -f 5)
BYTES_MEM=$((BYTES_MEM*4))
MY_TEMP_DIR=${TMPDIR:-/tmp}/CHC0009
mkdir $MY_TEMP_DIR

# echo Using $MY_TEMP_DIR for fifos and temporary bam index

rm -f $MY_TEMP_DIR/{a,b}.bam $MY_TEMP_DIR/{a,b,c}.sam
mkfifo $MY_TEMP_DIR/a.bam $MY_TEMP_DIR/b.bam $MY_TEMP_DIR/a.sam $MY_TEMP_DIR/b.sam $MY_TEMP_DIR/c.sam

/mnt/wigclust5/data/safe/kendall/bowtie-0.12.7/bowtie -S -t -n 2 -e 70 --chunkmbs 256 -3 0 -5 38 -m 1 --best --strata  hg19 <(gunzip -c /mnt/wigclust5/data/safe/kendall/riggs01/barcode.split/FC64VNGAAXX/s_5_sequence.WB1.txt.gz) 2> /mnt/wigclust5/data/safe/kendall/cox01/hg19.50k/mapped/CHC0009.bowtie.report.txt |
samtools view -Sbu -o - - |
(cd $MY_TEMP_DIR ; samtools sort -m ${BYTES_MEM} -o - -) |
samtools rmdup -s - - | tee $MY_TEMP_DIR/a.bam $MY_TEMP_DIR/b.bam > /mnt/wigclust5/data/safe/kendall/cox01/hg19.50k/mapped/CHC0009.rmdup.bam &
samtools index $MY_TEMP_DIR/a.bam &
samtools view -o - $MY_TEMP_DIR/b.bam | tee $MY_TEMP_DIR/a.sam $MY_TEMP_DIR/b.sam > $MY_TEMP_DIR/c.sam &
python /mnt/wigclust5/data/safe/kendall/cox01/hg19.50k/programs/varbin.50k.sam.py $MY_TEMP_DIR/a.sam /mnt/wigclust5/data/safe/kendall/cox01/hg19.50k/processed/CHC0009.varbin.50k.txt /mnt/wigclust5/data/safe/kendall/cox01/hg19.50k/processed/CHC0009.varbin.50k.stats.txt &
python /mnt/wigclust5/data/safe/kendall/cox01/hg19.50k/programs/varbin.5k.sam.py $MY_TEMP_DIR/b.sam /mnt/wigclust5/data/safe/kendall/cox01/hg19.50k/processed/CHC0009.varbin.5k.txt /mnt/wigclust5/data/safe/kendall/cox01/hg19.50k/processed/CHC0009.varbin.5k.stats.txt &
python /mnt/wigclust5/data/safe/kendall/cox01/hg19.50k/programs/varbin.20k.sam.py $MY_TEMP_DIR/c.sam /mnt/wigclust5/data/safe/kendall/cox01/hg19.50k/processed/CHC0009.varbin.20k.txt /mnt/wigclust5/data/safe/kendall/cox01/hg19.50k/processed/CHC0009.varbin.20k.stats.txt

echo waiting to finish
wait

rm $MY_TEMP_DIR/a.bam $MY_TEMP_DIR/b.bam $MY_TEMP_DIR/a.sam $MY_TEMP_DIR/b.sam $MY_TEMP_DIR/c.sam
mv $MY_TEMP_DIR/a.bam.bai /mnt/wigclust5/data/safe/kendall/cox01/hg19.50k/mapped/CHC0009.rmdup.bam.bai
rmdir $MY_TEMP_DIR/

exit
