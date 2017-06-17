library("DNAcopy", lib.loc="/mnt/wigclust5/data/safe/kendall/DNAcopy_1.26.0")

lowess.gc <- function(jtkx, jtky) {
        jtklow <- lowess(jtkx, log(jtky), f=0.05)
        jtkz <- approx(jtklow$x, jtklow$y, jtkx)
        return(exp(log(jtky) - jtkz$y))
}

cbs.segment01 <- function(indir, outdir, bad.bins, varbin.gc, varbin.data, sample.name, alt.sample.name, alpha, nperm, undo.SD, min.width) {
	gc <- read.table(varbin.gc, header=T)
	bad <- read.table(bad.bins, header=F)

	chrom.numeric <- substring(gc$bin.chrom, 4)
	chrom.numeric[which(gc$bin.chrom == "chrX")] <- "23"
	chrom.numeric[which(gc$bin.chrom == "chrY")] <- "24"
	chrom.numeric <- as.numeric(chrom.numeric)

	thisRatio <- read.table(paste(indir, varbin.data, sep="/"), header=F) 
	names(thisRatio) <- c("chrom", "chrompos", "abspos", "bincount", "ratio")
	thisRatio$chrom <- chrom.numeric
	a <- thisRatio$bincount + 1
	thisRatio$ratio <- a / mean(a)
	thisRatio$gc.content <- gc$gc.content
	thisRatio$lowratio <- lowess.gc(thisRatio$gc.content, thisRatio$ratio)

	a <- quantile(gc$bin.length, 0.985)
	thisRatioNobad <- thisRatio[which(bad[, 1] == 0), ]
	
	set.seed(25) 
	CNA.object <- CNA(log(thisRatio$lowratio, base=2), thisRatio$chrom, thisRatio$chrompos, data.type="logratio", sampleid=sample.name) 
	smoothed.CNA.object <- smooth.CNA(CNA.object) 
	segment.smoothed.CNA.object <- segment(smoothed.CNA.object, alpha=alpha, nperm=nperm, undo.splits="sdundo", undo.SD=undo.SD, min.width=min.width) 
	thisShort <- segment.smoothed.CNA.object[[2]]

	m <- matrix(data=0, nrow=nrow(thisRatio), ncol=1)	
	prevEnd <- 0
	for (i in 1:nrow(thisShort)) {
		thisStart <- prevEnd + 1
		thisEnd <- prevEnd + thisShort$num.mark[i]
		m[thisStart:thisEnd, 1] <- 2^thisShort$seg.mean[i]
		prevEnd = thisEnd
	}
	
	thisRatio$seg.mean.LOWESS <- m[, 1]

	chr <- thisRatio$chrom
	chr.shift <- c(chr[-1], chr[length(chr)])
	vlines <- c(1, thisRatio$abspos[which(chr != chr.shift) + 1], thisRatio$abspos[nrow(thisRatio)])
	hlines <- c(0.5, 1.0, 1.5, 2.0)
	chr.text <- c(1:22, "X", "Y")
	chr.text <- chr.text[c(1:14, 16, 18, 20, 22:24)]
	vlines.shift <- c(vlines[-1], 4*10^9)
	chr.at <- vlines + (vlines.shift - vlines) / 2
	chr.at <- chr.at[c(1:14, 16, 18, 20, 22:24)]
	x.at <- c(0, 0.5, 1, 1.5, 2, 2.5, 3) * 10^9
	x.labels <- c("0", "0.5", "1.0", "1.5", "2.0", "2.5", "3.0")
	y.at <- c(0.005, 0.020, 0.100, 0.500, 2.000)
	y.labels <- c("0.005", "0.020", "0.100", "0.500", "2.000")

	postscript(paste(outdir, "/", sample.name, ".wg.ps", sep=""), height=400, width=600)
	par(mar=c(5.1,4.1,4.1,4.1))
	plot(x=thisRatio$abspos, y=thisRatio$lowratio, log="y", main=paste(sample.name, alt.sample.name), xaxt="n", xlab="Genome Position Gb", yaxt="n", ylab="Ratio", col="#CCCCCC", cex=0.5)
	axis(1, at=x.at, labels=x.labels)
	axis(2, at=y.at, labels=y.labels)
	lines(x=thisRatio$abspos, y=thisRatio$lowratio, col="#CCCCCC", cex=0.5)
	points(x=thisRatio$abspos, y=thisRatio$seg.mean.LOWESS, col="#0000AA", cex=0.5)
	lines(x=thisRatio$abspos, y=thisRatio$seg.mean.LOWESS, col="#0000AA", cex=0.5)
	abline(h=hlines)
	abline(v=vlines)
	mtext(chr.text, at = chr.at)
	dev.off()
	
	write.table(thisRatio, sep="\t", file=paste(outdir, "/", sample.name, ".hg19.50k.k50.varbin.data.txt", sep=""), quote=F, row.names=F) 
	write.table(thisShort, sep="\t", file=paste(outdir, "/", sample.name, ".hg19.50k.k50.varbin.short.txt", sep=""), quote=F, row.names=F) 



	set.seed(25) 
	CNA.object <- CNA(log(thisRatioNobad$lowratio, base=2), thisRatioNobad$chrom, thisRatioNobad$chrompos, data.type="logratio", sampleid=sample.name) 
	smoothed.CNA.object <- smooth.CNA(CNA.object) 
	segment.smoothed.CNA.object <- segment(smoothed.CNA.object, alpha=alpha, nperm=nperm, undo.splits="sdundo", undo.SD=undo.SD, min.width=min.width) 
	thisShort <- segment.smoothed.CNA.object[[2]]

	m <- matrix(data=0, nrow=nrow(thisRatioNobad), ncol=1)	
	prevEnd <- 0
	for (i in 1:nrow(thisShort)) {
		thisStart <- prevEnd + 1
		thisEnd <- prevEnd + thisShort$num.mark[i]
		m[thisStart:thisEnd, 1] <- 2^thisShort$seg.mean[i]
		prevEnd = thisEnd
	}
	
	thisRatioNobad$seg.mean.LOWESS <- m[, 1]

	chr <- thisRatioNobad$chrom
	chr.shift <- c(chr[-1], chr[length(chr)])
	vlines <- c(1, thisRatioNobad$abspos[which(chr != chr.shift) + 1], thisRatioNobad$abspos[nrow(thisRatioNobad)])
	hlines <- c(0.5, 1.0, 1.5, 2.0)
	chr.text <- c(1:22, "X", "Y")
	chr.text <- chr.text[c(1:14, 16, 18, 20, 22:24)]
	vlines.shift <- c(vlines[-1], 4*10^9)
	chr.at <- vlines + (vlines.shift - vlines) / 2
	chr.at <- chr.at[c(1:14, 16, 18, 20, 22:24)]
	x.at <- c(0, 0.5, 1, 1.5, 2, 2.5, 3) * 10^9
	x.labels <- c("0", "0.5", "1.0", "1.5", "2.0", "2.5", "3.0")
	y.at <- c(0.005, 0.020, 0.100, 0.500, 2.000)
	y.labels <- c("0.005", "0.020", "0.100", "0.500", "2.000")

	postscript(paste(outdir, "/", sample.name, ".wg.nobad.ps", sep=""), height=400, width=600)
	par(mar=c(5.1,4.1,4.1,4.1))
	plot(x=thisRatioNobad$abspos, y=thisRatioNobad$lowratio, log="y", main=paste(sample.name, alt.sample.name), xaxt="n", xlab="Genome Position Gb", yaxt="n", ylab="Ratio", col="#CCCCCC", cex=0.5)
	axis(1, at=x.at, labels=x.labels)
	axis(2, at=y.at, labels=y.labels)
	lines(x=thisRatioNobad$abspos, y=thisRatioNobad$lowratio, col="#CCCCCC", cex=0.5)
	points(x=thisRatioNobad$abspos, y=thisRatioNobad$seg.mean.LOWESS, col="#0000AA", cex=0.5)
	lines(x=thisRatioNobad$abspos, y=thisRatioNobad$seg.mean.LOWESS, col="#0000AA", cex=0.5)
	abline(h=hlines)
	abline(v=vlines)
	mtext(chr.text, at = chr.at)
	dev.off()
	
	write.table(thisRatioNobad, sep="\t", file=paste(outdir, "/", sample.name, ".hg19.50k.k50.nobad.varbin.data.txt", sep=""), quote=F, row.names=F) 
	write.table(thisShort, sep="\t", file=paste(outdir, "/", sample.name, ".hg19.50k.k50.nobad.varbin.short.txt", sep=""), quote=F, row.names=F) 

}
cbs.segment01(indir=".", outdir=".", bad.bins="hg19.50k.k50.bad.bins.txt", varbin.gc="hg19.varbin.gc.content.50k.bowtie.k50.txt", varbin.data="SRR054616.varbin.50k.txt", sample.name="SRR054616", alt.sample.name="", alpha=0.05, nperm=1000, undo.SD=1.0, min.width=5)

