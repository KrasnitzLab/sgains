

## Install tools

```
conda config --add channels bioconda
conda install bowtie2 samtools bcftools
conda install wgsim

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


