all:


run: install
	python ../scpipe/r_pipeline.py

install:
	tar zcvf SCclust_0.1.7.tar.gz SCclust
	Rscript -e 'install.packages("SCclust_0.1.7.tar.gz", repos = NULL, type="source")'
