
ILLUMINA_ADAPTERS="/home/lubo/Local/anaconda3/envs/bioconda3/share/trimmomatic-0.36-3/adapters/TruSeq2-SE.fa"


CELLS="CJA0754
CJA0769
CJA0901
CJA0918
CJA1164
CJA1210
CJA1243
CJA1247
CJA1355"


for cell in $CELLS
do
  echo "Processing $cell"
  echo ">>> Trimming $cell"
  infile="$cell.fastq.gz"
  trimfile="$cell.jude58.fastq.gz"
  bamfile="$cell.jude58.bam"
  rmdupfile="$cell.jude58.rmdup.bam"

  gunzip -c $infile | \
    python adapter.clip02.min58.py GATCGGAAGAGCGG | \
    gzip > $trimfile
  echo ">>> Aligning $cell"
  bowtie -p 4 -S -t -n 2 -e 70 -3 18 -5 8 -m 1 --best --strata --solexa-quals \
    ../hg19/genomeindex <(gunzip -c $trimfile | head -n 10000) | \
    samtools view -Sbu -o $bamfile -
  echo ">>> Rmdup $cell"
  samtools sort -o - $bamfile | \
    samtools rmdup - $rmdupfile
  samtools index $rmdupfile
done

# trimmomatic SE CJA0918.fastq.gz CJA0918_trimmed.fastq.gz \
#    ILLUMINACLIP:$ILLUMINA_ADAPTERS:2:15:10
# 
# gunzip -c CJA0918.fastq.gz | python adapter.clip02.py GATCGGAAGAGCGG | gzip > CJA0918_jude02.fastq.gz
# 
# gunzip -c CJA0918.fastq.gz | python adapter.clip02.min58.py GATCGGAAGAGCGG | gzip > CJA0918_jude58.fastq.gz


