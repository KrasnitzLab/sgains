args<-commandArgs(trailingOnly = TRUE)
print(args)

library(parallel)
library(foreach)
library(doParallel)
cl<-parallel::makeCluster(2)
registerDoParallel(cl)
hgall<-c("chr1", "chr2")

# "chr3", "chr4", "chr5", "chr6", 
#	"chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
#	"chr13", "chr14", "chr15", 
#	"chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22",
#	"chrX", "chrY")


chrvec<-hgall

chrom <- "chr1";
command <- paste("python generate.reads.py ",args[2]," ",chrom,
		"|bowtie -S -t -v 0 -m 1 -f ",args[1],
		" - |python mappable.regions.py > data/readsim/",args[2],
		".txt 2> readsim/",args[2],".stderr",sep="");
print(command);


bychrom<-foreach(chrom=chrvec) %dopar% {
	command <- paste("python generate.reads.py ",args[2]," ",chrom,
			"|bowtie -S -t -v 0 -m 1 -f ",args[1],
			" - |python mappable.regions.py > data/readsim/",args[2],
			".txt 2> readsim/",args[2],".stderr",sep="");
	print(command);
	# system(command)
}
