args<-commandArgs(trailingOnly = TRUE)
library(parallel)
library(foreach)
library(doParallel)
mapfrom<-system(paste("ls raw/*",arg[1],sep=""),intern=T)
mapto<-substring(mapfrom,first=1,last=nchar(mapfrom)-nchar(arg[1]))
mapall<-foreach(mapme=mapfrom) %do% {
	system(paste("bowtie -S -t -n 2 -e 70 --chunkmbs 256 -3 18 -5 8 -m 1 --best --strata --solexa-quals genomeindex <(gunzip -c workhere/raw/",mapme," | python adapter.clip02.py GATCGGAAGAGCGG) | samtools view -Sbu -o workhere/mapped/",mapme,"bam -",sep=""))
}
