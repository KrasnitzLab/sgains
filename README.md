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

### Configure *s-GAINS* environment

In the root directory of the project there is a `setenv.sh` script. The purpose
of this script is to setup working environment of *s-GAINS*:

```

export PATH=$HOME/Local/anaconda3/bin:$PATH
source activate sgains3

export PATH=$(pwd)/tools:$PATH
export PYTHONPATH=$(pwd)/scpipe:$PYTHONPATH


```

The first line add *Anaconda 3* `bin` directory to the `PATH` variable,
so when you are using any program the *Anacond 3* `bin` directory whould be the
first directory to look for the tool. You may need to edit this line to point
to your local installation of *Anaconda 3*.

The second line activates prevously created anaconda environment.

The last two lines setup paths so that *s-GAINS* tools be accessible in your 
environment.

## Configure the *s-GAINS* pipeline

Example *s-GAINS* pipeline configuration:

```
genome:
    version: hg19
    pristine: data/hg19_safe
    work_dir: data/hg19
    index: genomeindex

reads:
    length: 100
    work_dir: data/R100
    bowtie_opts: ""
  
bins:
    bins_count: 10000
    work_dir: data/R100_B10k
    bins_boundaries: bins_boundaries.tsv

mapping:
    data_dir: data/test_study/raw
    data_glob: "*.fastq.gz"
    work_dir: data/test_study/bam
    bowtie_opts: "-S -t -n 2 -e 70 -3 18 -5 8 --solexa-quals"
    

varbin:
    data_dir: data/test_study/bam
    data_glob: "*.rmdup.bam"
    work_dir: data/test_study/try2
    suffix: "varbin.10k.txt"

segment:
    data_dir: data/test_study/try2
    data_glob: "*.varbin.10k.txt"
    work_dir: data/test_study/results
    study_name: "test_study"
```

Each section of this file configures different parts of ths *s-GAINS* pipeline.


## Usage of `sgains.py` tool

To interact with *s-GAINS* pipeline you invoke `sgains.py` command with different
parameters and subcommands. You can list available options of `sgains.py` using
`-h` option:

```
sgains.py -h
usage: sgains.py [-h] [-v] [-c path] [-n] [--force]
                 {mapping,varbin,segment} ...

sgains - sparse genomic analysis of individual nuclei by sequencing pipeline

USAGE

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         set verbosity level [default: 0]
  -c path, --config path
                        configuration file
  -n, --dry-run         perform a trial run with no changes made
  --force               allows overwriting nonempty results directory

subcommands:
  {mapping,varbin,segment}
    mapping             performs actual mapping of cell reads
    varbin              applies varbin algorithm to count read mappings in
                        each bin
    segment             segments bin counts and prepares the SCGV input data
```


