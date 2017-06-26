

## Install tools

```
conda config --add channels bioconda
conda install bowtie2 samtools bcftools
conda install wgsim

```

## Nature protocols paper

```
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5069701/
```

## Sam tools tutorial

```
http://biobits.org/samtools_primer.html
```


## Build the HG19 bowtie index

## Generate reads

```
cd data/hg19
python generate.reads.py 100 chr1 > ../readsim/chr01.reads.fa 
```

The script `generate.reads.py` uses about 10G RAM. 


## Align reads with `bowtie`

```
cd data/hg19
bowtie -S -t -v 0 -m 1 -f genomeindex ../readsim/chr01.reads.fa > ../readsim/chr01.reads.mapping.sam
```

## Mappable regions

```
python mappable.regions.py < data/readsim/chr01.reads.mapping.sam > data/readsim/chr01.stdout.txt 2> data/readsim/chr01.stderr.txt
```

## All mappable regions

```

```

## Sort mappable regions

```
sort -k 1,1 -k 2,2n data/readsim/chr01.stdout.txt > mappable.regions.sorted.txt
```


## Experiment with my tool chain

```
generate_mappable_regions.py -c scpipe.yml -l 100 -C chr1 > data/chr1.stdout.txt
```

```
generate_reads.py -c /mnt1/lubo/Work/single_cell_pipeline/scpipe.yml -C chr1 -l 100 | \
    bowtie2 -t -x genomeindex -f - | \
    mappable_regions.py -c /mnt1/lubo/Work/single_cell_pipeline/scpipe.yml
```

The original tool chain whould be something like this:

```
generate_reads_jc.py 100 hg19 | bowtie2 -t -x genomeindex -f - | mappable_regions_jc.py > chr1.orig.stdout.txt
```

## Time consumed to generate mappable regions

```
date && generate_mappable_regions.py -l 50 -c scpipe.yml -o data/readsim/all.50.mappable.txt && date
Fri Jun 16 14:14:26 EEST 2017
...
...
Sat Jun 17 04:20:49 EEST 2017
```

## Difference in results when writing reads

When we write reads with *biopython* `SeqIO` module  and when we write reads
by hand, the results of mappings is different. Not sure how bowtie works.



## FASTQ quality scores

```
https://en.wikipedia.org/wiki/FASTQ_format
```


## Illumina adapter trimming

* Trimmomatic: a flexible trimmer for Illumina sequence data
```
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4103590/
```

* Trimmomatic manual:
```
http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf
```

* Question: How To Best Deal With Adapter Contamination (Illumina)?
```
https://www.biostars.org/p/3461/
```

* Adapter and quality trimming of illumina data
```
http://www.ark-genomics.org/events-online-training-eu-training-course/adapter-and-quality-trimming-illumina-data
```

