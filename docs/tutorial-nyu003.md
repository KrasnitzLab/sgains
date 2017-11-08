# Example usage of `sGAINS` pipeline

## Download human reference genome HG19

Create a `data` directory where you plan to put your data files. Inside `data`
directory create a directory that would contain human reference genome `hg19_pristine`
and download a copy of the HG19 inside this directory. If you have `wget` installed
you can use following commands:

```
cd data/hg19_pristine
wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
tar zxvf chromFa.tar.gz
cd ..
```

## Build genome index

Go inside `data` directory and run following command to build the genome index

```
cd data
sgains.py genomeindex --genome-pristine hg19_pristine --genome-dir hg19
```

This command will create a subdirectory `data/hg19` that contains a copy of
HG19 reference genome with modified `chrX` with masked pseudo autosomal regions.
The `genomeindex` command will run `bowtie-build` to build genome index file 
that will be named `genomeindex.*`:

```
hg19/
├── genomeindex.1.ebwt
├── genomeindex.2.ebwt
├── genomeindex.3.ebwt
├── genomeindex.4.ebwt
├── genomeindex.rev.1.ebwt
└── genomeindex.rev.2.ebwt
```

## Preparation of mappable positions


### Precomputed mappable positions files

