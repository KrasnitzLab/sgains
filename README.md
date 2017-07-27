# Sparse Genomic Analysis of Individual Nuclei by Sequencing (s-GAINS)

## Anaconda environment setup

### Install Anacoda

### Install Anaconda 
* Go to anaconda web site 
[https://www.continuum.io/downloads](https://www.continuum.io/downloads)
and download the latest anaconda installer for your operating system. 

*s-GAINS* supports *Python 3.6* so you need to choos appropriate installer.
Note also that since *s-GAINS* uses *bioconda* channel the supported 
operating systems are only those supported for *biconda* (Linux and Mac OS X only
at the time of this writing).

* Install anaconda into suitable place on your local machine following
instructions from 
[https://docs.continuum.io/anaconda/install](https://docs.continuum.io/anaconda/install)


### Create `sgains` Anaconda environment

* After installing and activating *Anaconda* you need to create an environment to
use with `sgains` pipeline. To this end you need to use:

    ```
    conda create -n sgains
    source activate sgains
    ```

* After creating `sgains` environment you need to add *bioconda* and *r* channels:

    ```
    conda config --add channels bioconda
    conda config --add channels r
    ```

* Now you have to insall additional packages required by `sgains` tool:
    ```
    conda install bowtie
    conda install samtools bcftools trimmomatic biopython pysam fastqc
    conda install r-essentials
    conda install pandas numpy
    pip install python-box termcolor PyYAML pytest pytest-asyncio
    ```
