# Sparse Genomic Analysis of Individual Nuclei by Sequencing (s-GAINS)

**Please note that this software is still under construction and
not ready for general use**

## Anaconda environment setup

### Install Anaconda 
* Go to anaconda web site 
[https://www.continuum.io/downloads](https://www.continuum.io/downloads)
and download the latest anaconda installer for your operating system. 

* *s-GAINS* supports *Python 3.6* so you need to choose appropriate installer.
Note also that since *s-GAINS* uses *bioconda* channel the supported 
operating systems are only those supported for *bioconda* (at the time of
this writing that are Linux and Mac OS X).

* Install anaconda into suitable place on your local machine following
instructions from 
[https://docs.continuum.io/anaconda/install](https://docs.continuum.io/anaconda/install)


### Create `sgains` Anaconda environment

* After installing and activating *Anaconda* you need to create an environment to
use with `sgains` pipeline. To this end you need to use:

    ```
    conda create -n sgains3
    source activate sgains3
    ```

* After creating `sgains3` environment you need to add *bioconda* and *r* channels:

    ```
    conda config --add channels bioconda
    conda config --add channels r
    ```

* Now you have to insall additional packages required by `sgains.py` tool:

    ```
    conda install bowtie
    conda install samtools bcftools trimmomatic biopython pysam fastqc
    conda install r-essentials
    conda install pandas numpy
    pip install python-box termcolor PyYAML pytest pytest-asyncio
    ```

### Setup R environment

Go to `scripts` directory and invoke `setup.R` script:

```
cd scripts/
Rscript setup.R
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
so when you are using any program the *Anaconda's* `bin` directory whould be the
first directory to look for the tool. You may need to edit this line to point
to your local installation of *Anaconda 3*.

The second line activates prevously created anaconda environment.

The last two lines setup paths so that *s-GAINS* tools be accessible in your 
environment.


## Usage of `sgains.py` tool

To interact with *s-GAINS* pipeline you invoke `sgains.py` command with different
parameters and subcommands. You can list available options of `sgains.py` using
`-h` option:

```
sgains.py -h
usage: sgains.py [-h] [-v] [-c path] [-n] [--force]
                 {process,prepare,genomeindex,mappable-regions,bins,mapping,varbin,segment}
                 ...

sgains - sparse genomic analysis of individual nuclei by sequencing pipeline

USAGE

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         set verbosity level [default: 0]
  -c path, --config path
                        configuration file
  -n, --dry-run         perform a trial run with no changes made
  --force, -F           allows overwriting nonempty results directory

subcommands:
  {process,prepare,genomeindex,mappable-regions,bins,mapping,varbin,segment}
    process             combines mapping, varbin and segment subcommands into
                        single command
    prepare             combines all preparation steps (genomeindex,
                        mappable_regions, bins) into single command
    genomeindex         builds appropriate bowtie index for the reference genome
    mappable-regions    finds all mappable regions in specified genome
    bins                calculates all bins boundaries for specified bins
                        count and read length
    mapping             performs mapping of cell reads to reference genome
    varbin              applies varbin algorithm to count read mappings in
                        each bin
    segment             segments bin counts and prepares the SCGV input data
```

The `sgains.py` tool supports a list of common options:

* `--dry-run`, `-n` - this options instructs `sgains.py` to perform a trail run 
displaying information of commands that should be performed but without actualy 
running this commands

* `--force` - when `sgains.py` tool is run it checks if the result files or 
directories exists and if they already exists `sgains.py` stops whitout
performing any changes. To overwrite this behaivor you can use `--force` option

* `--config`, `-c` - instructs `sgains.py` which configuration file to use.

### Use of `process` subcommand

---

**Please note, that to use `process` subcommands 
(`mapping`, `varbin`, `segment` and `process`)
you need go through preparation steps.**

---

To list the options available for `process` subcommand use:

```
sgains.py process -h
usage: sgains.py process [-h] [--data-dir DATA_DIR] [--glob DATA_GLOB]
                         [--work-dir WORK_DIR] [--study-name STUDY_NAME]
                         [--genome-index GENOME_INDEX]
                         [--genome-dir GENOME_DIR]
                         [--genome-version GENOME_VERSION]
                         [--bowtie-opts BOWTIE_OPTS]
                         [--bins-boundaries BINS_BOUNDARIES]
                         [--bins-dir BINS_DIR]

optional arguments:
  -h, --help            show this help message and exit

input data options:
  --data-dir DATA_DIR, -i DATA_DIR
                        input data directory where the input data is located
                        (default: data/test_study/raw)
  --glob DATA_GLOB, -g DATA_GLOB
                        glob pattern for finding input data (default:
                        *.fastq.gz)

output data options:
  --work-dir WORK_DIR, -o WORK_DIR
                        output directory where results from processing are
                        stored (default: data/test_study/segment)
  --study-name STUDY_NAME, -s STUDY_NAME
                        study name (default: test_study)

genome index options:
  --genome-index GENOME_INDEX, -G GENOME_INDEX
                        genome index name (default: genomeindex)
  --genome-dir GENOME_DIR
                        genome index directory (default: data/hg19)
  --genome-version GENOME_VERSION
                        version of reference genome in use (supports only
                        hg19) (default: hg19)
  --bowtie-opts BOWTIE_OPTS
                        additional bowtie options (default: -p 10 -S -t -n 2
                        -e 70 -3 18 -5 8 --solexa-quals)

bins boundaries:
  --bins-boundaries BINS_BOUNDARIES, -B BINS_BOUNDARIES
                        bins boundaries filename (default:
                        hg19_R100_B10k_bins_boundaries.txt)
  --bins-dir BINS_DIR   bins working directory (default: data/R100_B10k)
```

* The data created by the `process` subcommand is grouped into a subdirectory, 
that is named with `--study-name` option. This name will be used when creating
results directories and when creating some of results files.

* The input for `process` subcommand are *FASTQ* files containing the reads for
each individual cell. All *FASTQ* files for given study are expected to be located
into single directory. You should specify this directory using `--data-dir` option.

* The results from `process` subcommand are stored into output data directory
which you should specify using `--work-dir` option. The process subcommand will
create directory named after the study name (passed with the `--study-name` option)
and inside that directory will create three additional subdirectories - `mapping`,
`varbin` and `segment`, that will contain intermediate results from respective
pipeline stages.

* The first `mapping` stage of the pipeline invokes `bowtie` to map reads from
*FASTQ* files. This stage needs a name of the bowtie index (user 
`--genome-index` option to specify bowtie index name) and a directory, 
where this index is located (use `--genome-dir` to pass this parameter). 

* When running `bowtie` for mapping *FASTQ* reads the pipeline passes a combination
of `-m 1 --best --strata`. If you need to pass additional options to `bowtie`
you can use `--bowtie-opts` option.

* The `varbin` stage of the pipeline needs a bins boundaries file prepared in
advance. You can pass bins boundaries file using `--bins-boundaries` option.


#### Example usage of `process` subcommand

If your input files are located into `data/test_study/raw` directory, then you
can use following command:

```
sgains.py process --data-dir data/test_study/raw \
    --work-dir data/try10 --study-name test_study \
    --genome-index genomeindex --genome-dir data/hg19/\
    --bins-boundaries data/R100_10k/bins_boundaries.tsv \
    --bowtie-opts "-S -t -n 2 -e 70 -3 18 -5 8 --solexa-quals"
```

The result of this command will be stored into directory with following structure:

    ```
    data/try10/
    └── test_study
        ├── mappings
        │   ├── ....rmdup.bam
        │   ├── ...
        │   └── ...
        ├── varbin
        │   ├── ....varbin.10k.txt
        │   ├── ...
        │   └── ...
        └── segment
            └── test_study
               ├── test_study.cells.txt
               ├── test_study.lowratio.quantal.R.ratio.csv
               ├── test_study.seg.quantal.R.seg.txt
               ├── test_study.smear1bpFisherPcloneTracks.clone.txt
               ├── test_study.smear1bpFisherTreePyP.tree.txt
               ├── test_study.smear1bpPinMat.featuremat.txt
               └── test_study.smear1bpPins.features.txt
    ```


### Usage of `mapping` subcommand

To list the options available for `mapping` subcommand use:

```
sgains.py mapping -h
usage: sgains.py mapping [-h] [--genome-index GENOME_INDEX]
                         [--genome-dir GENOME_DIR]
                         [--genome-version GENOME_VERSION]
                         [--data-dir DATA_DIR] [--glob DATA_GLOB]
                         [--bowtie-opts BOWTIE_OPTS] [--work-dir WORK_DIR]

optional arguments:
  -h, --help            show this help message and exit

genome index options:
  --genome-index GENOME_INDEX, -G GENOME_INDEX
                        genome index name (default: genomeindex)
  --genome-dir GENOME_DIR
                        genome index directory (default: data/hg19)
  --genome-version GENOME_VERSION
                        version of reference genome in use (supports only
                        hg19) (default: hg19)

input data options:
  --data-dir DATA_DIR, -i DATA_DIR
                        input data directory where the input data is located
                        (default: data/test_study/raw)
  --glob DATA_GLOB, -g DATA_GLOB
                        glob pattern for finding input data (default:
                        *.fastq.gz)
  --bowtie-opts BOWTIE_OPTS
                        additional bowtie options (default: -p 10 -S -t -n 2
                        -e 70 -3 18 -5 8 --solexa-quals)

output data options:
  --work-dir WORK_DIR, -o WORK_DIR
                        output directory where results from processing are
                        stored (default: data/kendall/bam)
```

### Use of `varbin` subcommand

To list the options available for `varbin` subcommand use:

```
sgains.py varbin -h
usage: sgains.py varbin [-h] [--data-dir DATA_DIR] [--glob DATA_GLOB]
                        [--work-dir WORK_DIR] [--suffix SUFFIX]
                        [--bins-boundaries BINS_BOUNDARIES]
                        [--bins-dir BINS_DIR]

optional arguments:
  -h, --help            show this help message and exit

input data options:
  --data-dir DATA_DIR, -i DATA_DIR
                        input data directory where the input data is located
                        (default: data/test_study/bam)
  --glob DATA_GLOB, -g DATA_GLOB
                        glob pattern for finding input data (default:
                        *.rmdup.bam)

output data options:
  --work-dir WORK_DIR, -o WORK_DIR
                        output directory where results from processing are
                        stored (default: data/test_study/try2)
  --suffix SUFFIX, -s SUFFIX
                        suffix for output files (default: varbin.10k.txt)

bins boundaries:
  --bins-boundaries BINS_BOUNDARIES, -B BINS_BOUNDARIES
                        bins boundaries filename (default:
                        hg19_R100_B10k_bins_boundaries.txt)
  --bins-dir BINS_DIR   bins working directory (default: data/R100_B10k)
```



### Use of `segment` subcommand

To list options available for `segment` subcommand use:

```
sgains.py segment -h
usage: sgains.py segment [-h] [--data-dir DATA_DIR] [--glob DATA_GLOB]
                         [--work-dir WORK_DIR] [--study-name STUDY_NAME]
                         [--bins-boundaries BINS_BOUNDARIES]
                         [--bins-dir BINS_DIR]

optional arguments:
  -h, --help            show this help message and exit

input data options:
  --data-dir DATA_DIR, -i DATA_DIR
                        input data directory where the input data is located
                        (default: data/test_study/try2)
  --glob DATA_GLOB, -g DATA_GLOB
                        glob pattern for finding input data (default:
                        *.varbin.10k.txt)

output data options:
  --work-dir WORK_DIR, -o WORK_DIR
                        output directory where results from processing are
                        stored (default: data/test_study/segment)
  --study-name STUDY_NAME, -s STUDY_NAME
                        study name (default: test_study)

bins boundaries:
  --bins-boundaries BINS_BOUNDARIES, -B BINS_BOUNDARIES
                        bins boundaries filename (default:
                        hg19_R100_B10k_bins_boundaries.txt)
  --bins-dir BINS_DIR   bins working directory (default: data/R100_B10k)
```

## Pipeline preparation

### Usage of `genomeindex` subcommand
The `genomeindex` subcommand builds the bowtie index for the reference genome. To
list the available options use:

```
sgains.py genomeindex -h
usage: sgains.py genomeindex [-h] [--genome-index GENOME_INDEX]
                             [--genome-dir GENOME_DIR]
                             [--genome-version GENOME_VERSION]
                             [--genome-pristine GENOME_PRISTINE]

optional arguments:
  -h, --help            show this help message and exit

genome index options:
  --genome-index GENOME_INDEX, -G GENOME_INDEX
                        genome index name (default: genomeindex)
  --genome-dir GENOME_DIR
                        genome index directory (default: data/hg19)
  --genome-version GENOME_VERSION
                        version of reference genome in use (supports only
                        hg19) (default: hg19)
  --genome-pristine GENOME_PRISTINE
                        directory where clean copy of reference genome is
                        located (default: data/hg19_safe)
```

### Usage of `mappable-regions` subcommand

This command find all uniquely mappable regions of the reference genome with
given length.

This step is computationally expesive and could take days in CPU time. 

To save this step you can use files with precomputed mappable regions that
could be found at:

* For Human Reference Genome **HG19** with read length **100bp**: 
[hg19_R100_mappable_regions.txt.gz](https://github.com/KrasnitzLab/sgains/releases/download/1.0/hg19_R100_mappable_regions.txt.gz)


You can download and unzip some of these files and use them into following
stages of the pipeline preparation.

If you want to build your own mappable regions file you can use `mappable-regions`
subcomman. To run this command you will need genome index build from `genomeindex`
subommand.

To list the options available for this subcommand use:

```
sgains.py mappable-regions -h
usage: sgains.py mappable-regions [-h] [--mappable-dir MAPPABLE_DIR]
                                  [--mappable-regions MAPPABLE_REGIONS]
                                  [--read-length LENGTH]
                                  [--bowtie-opts BOWTIE_OPTS]
                                  [--genome-index GENOME_INDEX]
                                  [--genome-dir GENOME_DIR]
                                  [--genome-version GENOME_VERSION]

optional arguments:
  -h, --help            show this help message and exit

mappable regions options:
  --mappable-dir MAPPABLE_DIR, -m MAPPABLE_DIR
                        directory where mappable regions file is stroed
                        (default: data/R100)
  --mappable-regions MAPPABLE_REGIONS, -M MAPPABLE_REGIONS
                        filename where mappable regions are stored (default:
                        hg19_R100_mappable_regions.txt)
  --read-length LENGTH, -l LENGTH
                        read length to use for generation of mappable regions
                        (default: 100)
  --bowtie-opts BOWTIE_OPTS
                        additional bowtie options (default: -p 20)

genome index options:
  --genome-index GENOME_INDEX, -G GENOME_INDEX
                        genome index name (default: genomeindex)
  --genome-dir GENOME_DIR
                        genome index directory (default: data/hg19)
  --genome-version GENOME_VERSION
                        version of reference genome in use (supports only
                        hg19) (default: hg19)
```


### Usage of `bins` subcommand

The `bins` subcommand computes the bins boudaries.

To list options available for `bins` subcommand use:

```
sgains.py bins -h
usage: sgains.py bins [-h] [--mappable-dir MAPPABLE_DIR]
                      [--mappable-regions MAPPABLE_REGIONS]
                      [--bins-boundaries BINS_BOUNDARIES]
                      [--bins-dir BINS_DIR] [--bins-count BINS_COUNT]

optional arguments:
  -h, --help            show this help message and exit

mappable regions options:
  --mappable-dir MAPPABLE_DIR, -m MAPPABLE_DIR
                        directory where mappable regions file is stroed
                        (default: data/R100)
  --mappable-regions MAPPABLE_REGIONS, -M MAPPABLE_REGIONS
                        filename where mappable regions are stored (default:
                        hg19_R100_mappable_regions.txt)

bins boundaries:
  --bins-boundaries BINS_BOUNDARIES, -B BINS_BOUNDARIES
                        bins boundaries filename (default:
                        hg19_R100_B10k_bins_boundaries.txt)
  --bins-dir BINS_DIR   bins working directory (default: data/R100_B10k)
  --bins-count BINS_COUNT, -C BINS_COUNT
                        number of bins (default: 10000)
```

### Usage of `prepare` subcommand

This command combines all preparation subcommands into one. Please note, this
is very computationaly expensive command an use it with care.

To list options available for `prepare` subcommand use:

```
sgains.py prepare -h
usage: sgains.py prepare [-h] [--genome-index GENOME_INDEX]
                         [--genome-dir GENOME_DIR]
                         [--genome-version GENOME_VERSION]
                         [--genome-pristine GENOME_PRISTINE]
                         [--mappable-dir MAPPABLE_DIR]
                         [--mappable-regions MAPPABLE_REGIONS]
                         [--read-length LENGTH] [--bowtie-opts BOWTIE_OPTS]
                         [--bins-boundaries BINS_BOUNDARIES]
                         [--bins-dir BINS_DIR] [--bins-count BINS_COUNT]

optional arguments:
  -h, --help            show this help message and exit

genome index options:
  --genome-index GENOME_INDEX, -G GENOME_INDEX
                        genome index name (default: genomeindex)
  --genome-dir GENOME_DIR
                        genome index directory (default: data/hg19)
  --genome-version GENOME_VERSION
                        version of reference genome in use (supports only
                        hg19) (default: hg19)
  --genome-pristine GENOME_PRISTINE
                        directory where clean copy of reference genome is
                        located (default: data/hg19_safe)

mappable regions options:
  --mappable-dir MAPPABLE_DIR, -m MAPPABLE_DIR
                        directory where mappable regions file is stroed
                        (default: data/R100)
  --mappable-regions MAPPABLE_REGIONS, -M MAPPABLE_REGIONS
                        filename where mappable regions are stored (default:
                        hg19_R100_mappable_regions.txt)
  --read-length LENGTH, -l LENGTH
                        read length to use for generation of mappable regions
                        (default: 100)
  --bowtie-opts BOWTIE_OPTS
                        additional bowtie options (default: -p 20)

bins boundaries:
  --bins-boundaries BINS_BOUNDARIES, -B BINS_BOUNDARIES
                        bins boundaries filename (default:
                        hg19_R100_B10k_bins_boundaries.txt)
  --bins-dir BINS_DIR   bins working directory (default: data/R100_B10k)
  --bins-count BINS_COUNT, -C BINS_COUNT
                        number of bins (default: 10000)
```



## Configure the *s-GAINS* pipeline

Example *s-GAINS* pipeline configuration:

```
genome:
    version: hg19
    data_dir: data/hg19_safe
    work_dir: data/hg19
    genome_index: genomeindex

mappable_regions:
    length: 100
    work_dir: data/R100
    bowtie_opts: "-p 20"
  
bins:
    bins_count: 10000
    work_dir: data/R100_B10k
    bins_boundaries: hg19_R100_B10k_bins_boundaries.txt

mapping:
    data_dir: data/test_study/raw
    data_glob: "*.fastq.gz"
    work_dir: data/kendall/bam
    bowtie_opts: "-p 10 -S -t -n 2 -e 70 -3 18 -5 8 --solexa-quals"
    

varbin:
    data_dir: data/test_study/bam
    data_glob: "*.rmdup.bam"
    work_dir: data/test_study/try2
    suffix: "varbin.10k.txt"

segment:
    data_dir: data/test_study/try2
    data_glob: "*.varbin.10k.txt"
    work_dir: data/test_study/segment
    study_name: "test_study"
```

Each section of this file configures different parts of ths *s-GAINS* pipeline.
