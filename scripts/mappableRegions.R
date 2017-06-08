args<-commandArgs(trailingOnly = TRUE)
library(parallel)
library(foreach)
library(doParallel)
cl<-parallel::makeCluster(getOption("cl.cores",detectCores()))
registerDoParallel(cl)
hgall<-c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", 
	"chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", 
	"chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
if(arg[2]=="hg")chrvec<-hgall
else if(arg[2]=="hgdm")chrvec<-c(hgall,"dm_chrX","dm_chr2L","dm_chr2R",
	"dm_chr3L","dm_chr3R","dm_chr4")
bychrom<-foreach(chrom=chrvec) %dopar% {
	system(paste("python generate.reads.py ",arg[1]," ",chrom,
		"|bowtie -S -t -v 0 -m 1 -f ",arg[2],
		" - |python mappable.regions.py > readsim/",arg[2],
		".txt 2> readsim/",arg[2],".stderr",sep=""))
}
