all:


run: install
	python ../scpipe/r_pipeline.py

install:
	tar zcvf SCclust_0.1.5.tar.gz SCclust
	Rscript -e 'install.packages("SCclust_0.1.5.tar.gz", repos = NULL, type="source")'
