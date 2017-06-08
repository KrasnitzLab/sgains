# Define base image
FROM ubuntu:latest

# Install required packages
RUN apt-get update && apt-get install -y \
 	  wget \
	autoconf \
	automake \
	make \
	g++ \
	gcc \
	build-essential \ 
	zlib1g-dev \
	libgsl0-dev \
	python \
	python-pip \
	curl \
	git \
	r-base \
	bowtie \
	samtools


# Change workdir
WORKDIR /opt
#install python-config (missing)
RUN pip install python-config
#install DNAcopy, TBEST and SCclust for R
RUN echo 'source("https://bioconductor.org/biocLite.R")' > Rtools.R
RUN echo 'biocLite("DNAcopy")' >> Rtools.R
RUN echo 'install.packages("TBEST",repos="http://cran.us.r-project.org")' >> Rtools.R
RUN echo 'install.packages("rPython",repos="http://cran.us.r-project.org")' >> Rtools.R
RUN R CMD BATCH Rtools.R
#copy and install the in-house package
COPY SCclust_0.1.4.tar.gz /opt/
RUN R CMD INSTALL SCclust_0.1.4.tar.gz

#copy Python, R and shell scripts
COPY hg19.chrY.psr.py /opt/
#COPY hg38.chrY.psr.py /opt/
COPY generate.reads.py /opt/
COPY bin.boundaries.py /opt/
COPY chrom.mappable.bowtie.py /opt/
COPY mappable.regions.py /opt/
COPY varbin.gc.content.bowtie.py /opt/
COPY adapter.clip02.py /opt/
COPY build.index.bash /opt/
COPY chrYmask.R /opt/
COPY mappableRegions.R /opt/
COPY serialMapper.R /opt/


COPY Dockerfile /opt/

# Maintainer
MAINTAINER Alex Krasnitz, Cold Spring Harbor Laboratory <krasnitz@cshl.edu>
