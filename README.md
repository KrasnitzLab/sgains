# Sparse Genomic Analysis of Individual Nuclei by Sequencing (s-GAINS)

[![DOI](https://zenodo.org/badge/98500084.svg)](https://zenodo.org/badge/latestdoi/98500084)

This document describes how to setup `s-GAINS` pipeline tool and its basic command.

Short tutorial on how to use this tool could be found in
[Example usage of `sGAINS` pipeline](docs/tutorial-navin2011.md)

## Anaconda environment setup

### Install Anaconda

* Go to anaconda web site
    [https://www.anaconda.com/distribution/](https://www.anaconda.com/distribution/)
    and download the latest anaconda installer for your operating system.

* *s-GAINS* supports *Python 3.7* or greater so you need to choose an
    appropriate installer. Note also that since *s-GAINS* uses *bioconda*
    channel the supported operating systems are only those supported for
    *bioconda* (at the time of this writing these are Linux and Mac OS X).

* Install anaconda into suitable place on your local machine following
    instructions from
    [https://docs.anaconda.com/anaconda/install/](https://docs.anaconda.com/anaconda/install/)

### Create `sgains` Anaconda environment

* After installing and activating *Anaconda* you need to create an environment to
    use with `sgains` pipeline. To this end you need to use:

    ```bash
    conda create -n sgains3
    conda activate sgains3
    ```

### Install `sgains` anaconda package

* *sGAINS* tools are distributed as a conda package through `krasnitzlab`
    Annaconda channel. So to install *sGAINS* tools use:

    ```bash
    conda install -c defaults -c conda-forge -c krasnitzlab -c bioconda sgains
    ```

    This command should install all the packages and tools need for
    proper functioning of `sgains-tools`.

* After this command finishes, you should be able to use
    `sgains-tools` command:

    ```bash
    sgains-tools --help
    ```

### Install SCGV viewer package

To visualize results of `sgains-tools` you may need `SCGV` viewer.

`SCGV` package is available from `KrasnitzLab` Anaconda channel.
You can to install it using using following command:

```bash
conda install -c krasnitzlab scgv
```

## Usage of sgains docker container

Instead of seting up `sgains` environment you can use `krasnitzlab/sgains`
docker container image to run the pipeline. To this end you need to have *Docker*
tools installed and configured on your computer (please look for instructions
in the official [*Docker* documentation](https://docs.docker.com).

### Download *s-GAINS* container image

Once you have Docker installed and configured you can pull `krasnitzlab/sgains`
docker container image by using docker pull command:

```bash
docker pull krasnitzlab/sgains
```

### Run *s-GAINS* container in interactive mode

You can run the `sgains` container interactively by using:

```bash
docker run -i -v /data/pathname:/data -t krasnitzlab/sgains /bin/bash
```

where `/data/pathname` is a full pathname to a folder on your local machine,
where data you want to process is located.

### Run *s-GAINS* commands

You can use this docker container to run all subcommans of
`sgains-tools` using
following sintax:

```bash
docker run -i -v /data/pathname:/data -t krasnitzlab/sgains sgains-tools <arg1> <arg2> ...
```

In this way you can run any `sgains-tools` subcommand with appropriate arguments
you need.

## Usage of `sgains-tools` tool

To interact with *s-GAINS* pipeline you invoke `sgains-tools` command with different
parameters and subcommands. You can list available options of `sgains-tools` using
`-h` option:

```bash
sgains-tools -h
usage: sgains-tools [-h] [-v] [-c path] [-n] [--force] [--parallel PARALLEL]
                    [--sge]
                    {genome,mappable-regions,bins,prepare,mapping,extract-10x,varbin,varbin-10x,scclust,process}
                    ...

sgains - sparse genomic analysis of individual nuclei by sequencing pipeline
USAGE

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         set verbosity level [default: 0]
  -c path, --config path
                        configuration file (default: None)
  -n, --dry-run         perform a trial run with no changes made (default:
                        False)
  --force, -F           allows overwriting nonempty results directory
                        (default: False)
  --parallel PARALLEL, -p PARALLEL
                        number of task to run in parallel (default: 1)
  --sge                 parallelilizes commands using SGE cluster manager
                        (default: False)

sGAINS subcommands:
  {genome,mappable-regions,bins,prepare,mapping,extract-10x,varbin,varbin-10x,scclust,process}
    genomeindex         builds appropriate hisat2 or bowtie index for the
                        reference genome
    mappable-regions    finds all mappable regions in specified genome
    bins                calculates all bins boundaries for specified bins
                        count and read length
    prepare             combines all preparation steps ('genome', 'mappable-
                        regions' and 'bins') into single command
    mapping             performs mapping of cells reads to the reference
                        genome
    extract-10x         extracts cells reads from 10x Genomics datasets
    varbin              applies varbin algorithm to count read mappings in
                        each bin
    varbin-10x          applies varbin algorithm to count read mappings in
                        each bin to 10x Genomics datasets without realigning
    scclust             segmentation and clustering based bin counts and
                        preparation of the SCGV input data
    process             combines all process steps ('mapping', 'varbin' and
                        'scclust') into single command
```

The `sgains-tools` tool supports a list of common options:

* `--dry-run`, `-n` - this option instructs `sgains-tools` to perform a trail run
    displaying information of commands that should be performed but without actualy
    running these commands

* `--force` - when `sgains-tools` tool is run it checks if the result files or
    directories already exist and, if they do, `sgains-tools` stops whitout
    making any changes. To override this behaivor you can use the  `--force` option

* `--config`, `-c` - instructs `sgains-tools` which configuration file to use.

* `--parallel`, `-p` - instructs `sgains-tools` to parallelize work on subcommands
    called.

* `--sge` - parallellilize execution using SGE.


## Pipeline preparation

### Usage of `genomeindex` subcommand

The `genomeindex` subcommand builds the bowtie index for the reference genome. To
list the available options use:

```bash
sgains-tools genomeindex -h
usage: sgains-tools genomeindex [-h] [--aligner-name ALIGNER_NAME]
                                [--genome-version GENOME_VERSION]
                                [--genome-pristine-dir GENOME_PRISTINE_DIR]
                                [--genome-dir GENOME_DIR]
                                [--genomeindex-prefix GENOMEINDEX_PREFIX]

optional arguments:
  -h, --help            show this help message and exit

aligner group::
  --aligner-name ALIGNER_NAME
                        aligner to use in sGAINS subcommands (default: bowtie)

genome group::
  --genome-version GENOME_VERSION
                        version of reference genome to use (default: hg19)
  --genome-pristine-dir GENOME_PRISTINE_DIR
                        directory where clean copy of reference genome is
                        located (default: None)
  --genome-dir GENOME_DIR
                        genome index working directory (default: None)
  --genomeindex-prefix GENOMEINDEX_PREFIX
                        genome index prefix (default: genomeindex)
```

### Usage of `mappable-regions` subcommand

This command find all uniquely mappable regions of the reference genome with
given length.

This step is computationally expesive and could take days in CPU time.

To save this step you can use files with precomputed mappable regions that
could be found at:

* For Human Reference Genome **HG19** with read length **50bp**:
    [hg19_R50_mappable_regions.txt.gz](https://github.com/KrasnitzLab/sgains/releases/download/1.0.0RC1/hg19_R50_mappable_regions.txt.gz)

You can download and unzip some of these files and use them into following
stages of the pipeline preparation.

If you want to build your own mappable regions file you can use `mappable-regions`
subcommand. To run this command you will need genome index build from `genomeindex`
subommand.

To list the options available for this subcommand use:

```bash
sgains-tools mappable-regions -h
usage: sgains-tools mappable-regions [-h] [--aligner-name ALIGNER_NAME]
                                     [--genome-version GENOME_VERSION]
                                     [--genome-pristine-dir GENOME_PRISTINE_DIR]
                                     [--genome-dir GENOME_DIR]
                                     [--genomeindex-prefix GENOMEINDEX_PREFIX]
                                     [--mappable-read-length MAPPABLE_READ_LENGTH]
                                     [--mappable-dir MAPPABLE_DIR]
                                     [--mappable-file MAPPABLE_FILE]
                                     [--mappable-aligner-options MAPPABLE_ALIGNER_OPTIONS]

optional arguments:
  -h, --help            show this help message and exit

aligner group::
  --aligner-name ALIGNER_NAME
                        aligner to use in sGAINS subcommands (default: bowtie)

genome group::
  --genome-version GENOME_VERSION
                        version of reference genome to use (default: hg19)
  --genome-pristine-dir GENOME_PRISTINE_DIR
                        directory where clean copy of reference genome is
                        located (default: None)
  --genome-dir GENOME_DIR
                        genome index working directory (default: None)
  --genomeindex-prefix GENOMEINDEX_PREFIX
                        genome index prefix (default: genomeindex)

mappable_regions group::
  --mappable-read-length MAPPABLE_READ_LENGTH
                        read length to use for generation of mappable regions
                        (default: 100)
  --mappable-dir MAPPABLE_DIR
                        directory where mappable regions working files are
                        stored (default: None)
  --mappable-file MAPPABLE_FILE
                        filename for mappable regions results (default:
                        mappable_regions.txt)
  --mappable-aligner-options MAPPABLE_ALIGNER_OPTIONS
                        additional aligner options for use when computing
                        uniquely mappable regions (default: )
```

### Usage of `bins` subcommand

The `bins` subcommand computes the bins boudaries.

To list options available for `bins` subcommand use:

```bash
sgains-tools bins -h
usage: sgains-tools bins [-h] [--genome-version GENOME_VERSION]
                         [--genome-pristine-dir GENOME_PRISTINE_DIR]
                         [--genome-dir GENOME_DIR]
                         [--genomeindex-prefix GENOMEINDEX_PREFIX]
                         [--mappable-read-length MAPPABLE_READ_LENGTH]
                         [--mappable-dir MAPPABLE_DIR]
                         [--mappable-file MAPPABLE_FILE]
                         [--mappable-aligner-options MAPPABLE_ALIGNER_OPTIONS]
                         [--bins-count BINS_COUNT] [--bins-dir BINS_DIR]
                         [--bins-file BINS_FILE]

optional arguments:
  -h, --help            show this help message and exit

genome group::
  --genome-version GENOME_VERSION
                        version of reference genome to use (default: hg19)
  --genome-pristine-dir GENOME_PRISTINE_DIR
                        directory where clean copy of reference genome is
                        located (default: None)
  --genome-dir GENOME_DIR
                        genome index working directory (default: None)
  --genomeindex-prefix GENOMEINDEX_PREFIX
                        genome index prefix (default: genomeindex)

mappable_regions group::
  --mappable-read-length MAPPABLE_READ_LENGTH
                        read length to use for generation of mappable regions
                        (default: 100)
  --mappable-dir MAPPABLE_DIR
                        directory where mappable regions working files are
                        stored (default: None)
  --mappable-file MAPPABLE_FILE
                        filename for mappable regions results (default:
                        mappable_regions.txt)
  --mappable-aligner-options MAPPABLE_ALIGNER_OPTIONS
                        additional aligner options for use when computing
                        uniquely mappable regions (default: )

bins group::
  --bins-count BINS_COUNT
                        number of bins (default: 10000)
  --bins-dir BINS_DIR   bins working directory (default: None)
  --bins-file BINS_FILE
                        bins boundaries filename (default:
                        bins_boundaries.txt)
```

## Processing sequence data

### Use of `process` subcommand

---

**Please note, that to use `process` subcommands
(`mapping`, `varbin`, `scclust` and `process`)
you need go through all the preparation steps.**

---

To list the options available for `process` subcommand use:

```bash
sgains-tools process -h
usage: sgains-tools process [-h] [--aligner-name ALIGNER_NAME]
                            [--genome-version GENOME_VERSION]
                            [--genome-pristine-dir GENOME_PRISTINE_DIR]
                            [--genome-dir GENOME_DIR]
                            [--genomeindex-prefix GENOMEINDEX_PREFIX]
                            [--reads-dir READS_DIR]
                            [--reads-suffix READS_SUFFIX]
                            [--mapping-dir MAPPING_DIR]
                            [--mapping-suffix MAPPING_SUFFIX]
                            [--mapping-aligner-options MAPPING_ALIGNER_OPTIONS]
                            [--bins-count BINS_COUNT] [--bins-dir BINS_DIR]
                            [--bins-file BINS_FILE] [--varbin-dir VARBIN_DIR]
                            [--varbin-suffix VARBIN_SUFFIX]
                            [--scclust-case SCCLUST_CASE]
                            [--scclust-dir SCCLUST_DIR]
                            [--scclust-cytoband-file SCCLUST_CYTOBAND_FILE]
                            [--scclust-nsim SCCLUST_NSIM]
                            [--scclust-sharemin SCCLUST_SHAREMIN]
                            [--scclust-fdrthres SCCLUST_FDRTHRES]
                            [--scclust-nshare SCCLUST_NSHARE]
                            [--scclust-climbtoshare SCCLUST_CLIMBTOSHARE]

optional arguments:
  -h, --help            show this help message and exit

aligner group::
  --aligner-name ALIGNER_NAME
                        aligner to use in sGAINS subcommands (default: bowtie)

genome group::
  --genome-version GENOME_VERSION
                        version of reference genome to use (default: hg19)
  --genome-pristine-dir GENOME_PRISTINE_DIR
                        directory where clean copy of reference genome is
                        located (default: None)
  --genome-dir GENOME_DIR
                        genome index working directory (default: None)
  --genomeindex-prefix GENOMEINDEX_PREFIX
                        genome index prefix (default: genomeindex)

reads group::
  --reads-dir READS_DIR
                        data directory where sequencing reads are located
                        (default: None)
  --reads-suffix READS_SUFFIX
                        reads files suffix pattern (default: .fastq.gz)

mapping group::
  --mapping-dir MAPPING_DIR
                        data directory where mapping files are located
                        (default: None)
  --mapping-suffix MAPPING_SUFFIX
                        mapping files suffix pattern (default: .rmdup.bam)
  --mapping-aligner-options MAPPING_ALIGNER_OPTIONS
                        additional aligner mapping options (default: )

bins group::
  --bins-count BINS_COUNT
                        number of bins (default: 10000)
  --bins-dir BINS_DIR   bins working directory (default: None)
  --bins-file BINS_FILE
                        bins boundaries filename (default:
                        bins_boundaries.txt)

varbin group::
  --varbin-dir VARBIN_DIR
                        varbin working directory (default: None)
  --varbin-suffix VARBIN_SUFFIX
                        varbin files suffix pattern (default: .varbin.txt)

scclust group::
  --scclust-case SCCLUST_CASE
                        SCclust case name (default: None)
  --scclust-dir SCCLUST_DIR
                        SCclust working directory (default: None)
  --scclust-cytoband-file SCCLUST_CYTOBAND_FILE
                        location of cyto band description file (default: None)
  --scclust-nsim SCCLUST_NSIM
                        SCclust number of simulations (default: 150)
  --scclust-sharemin SCCLUST_SHAREMIN
                        SCclust sharemin parameter (default: 0.85)
  --scclust-fdrthres SCCLUST_FDRTHRES
                        SCclust fdrthres parameter (default: -3)
  --scclust-nshare SCCLUST_NSHARE
                        SCclust nshare parameter (default: 4)
  --scclust-climbtoshare SCCLUST_CLIMBTOSHARE
                        SCclust climbtoshare parameter (default: 5)
```

* The data created by the `process` subcommand are placed in a subdirectory,
    whose name is specified with `--output-dir` option. This name will be used
    when creating
    the result directory structure.

* The input for `process` subcommand are *FASTQ* files containing the reads for
    each individual cell. All *FASTQ* files for given study are expected to be
    located into single directory. You should specify this directory using
    `--reads-dir` option.

* The results from `process` subcommand are stored in the output data directory,
    as specified using `--output-dir` option. The process subcommand will
    create a directory and inside that directory it will create three additional
    subdirectories - `mapping`,
    `varbin` and `scclust`. These will contain intermediate results from the
    respective pipeline stages.

* The first `mapping` stage of the pipeline invokes `bowtie` to map reads from
    *FASTQ* files. This stage needs a name of the bowtie index (user
    `--genome-index` option to specify bowtie index name) and a directory,
    where this index is located (use `--genome-dir` to pass this parameter).

* If you need to pass additional options to `bowtie` to control mapping reads
    you can use `--mapping-bowtie-opts` option.

* The `varbin` stage of the pipeline needs a bins boundaries file prepared in
    advance. You can pass bins boundaries file using `--bins-boundaries` option.

### Usage of `mapping` subcommand

To list the options available for `mapping` subcommand use:

```bash
sgains-tools mapping -h
usage: sgains-tools mapping [-h] [--aligner-name ALIGNER_NAME]
                            [--genome-version GENOME_VERSION]
                            [--genome-pristine-dir GENOME_PRISTINE_DIR]
                            [--genome-dir GENOME_DIR]
                            [--genomeindex-prefix GENOMEINDEX_PREFIX]
                            [--reads-dir READS_DIR]
                            [--reads-suffix READS_SUFFIX]
                            [--mapping-dir MAPPING_DIR]
                            [--mapping-suffix MAPPING_SUFFIX]
                            [--mapping-aligner-options MAPPING_ALIGNER_OPTIONS]

optional arguments:
  -h, --help            show this help message and exit

aligner group::
  --aligner-name ALIGNER_NAME
                        aligner to use in sGAINS subcommands (default: bowtie)

genome group::
  --genome-version GENOME_VERSION
                        version of reference genome to use (default: hg19)
  --genome-pristine-dir GENOME_PRISTINE_DIR
                        directory where clean copy of reference genome is
                        located (default: None)
  --genome-dir GENOME_DIR
                        genome index working directory (default: None)
  --genomeindex-prefix GENOMEINDEX_PREFIX
                        genome index prefix (default: genomeindex)

reads group::
  --reads-dir READS_DIR
                        data directory where sequencing reads are located
                        (default: None)
  --reads-suffix READS_SUFFIX
                        reads files suffix pattern (default: .fastq.gz)

mapping group::
  --mapping-dir MAPPING_DIR
                        data directory where mapping files are located
                        (default: None)
  --mapping-suffix MAPPING_SUFFIX
                        mapping files suffix pattern (default: .rmdup.bam)
  --mapping-aligner-options MAPPING_ALIGNER_OPTIONS
                        additional aligner mapping options (default: )
```

### Use of `varbin` subcommand

To list the options available for `varbin` subcommand use:

```bash
sgains-tools varbin -h
usage: sgains-tools varbin [-h] [--bins-count BINS_COUNT]
                           [--bins-dir BINS_DIR] [--bins-file BINS_FILE]
                           [--mapping-dir MAPPING_DIR]
                           [--mapping-suffix MAPPING_SUFFIX]
                           [--mapping-aligner-options MAPPING_ALIGNER_OPTIONS]
                           [--varbin-dir VARBIN_DIR]
                           [--varbin-suffix VARBIN_SUFFIX]

optional arguments:
  -h, --help            show this help message and exit

bins group::
  --bins-count BINS_COUNT
                        number of bins (default: 10000)
  --bins-dir BINS_DIR   bins working directory (default: None)
  --bins-file BINS_FILE
                        bins boundaries filename (default:
                        bins_boundaries.txt)

mapping group::
  --mapping-dir MAPPING_DIR
                        data directory where mapping files are located
                        (default: None)
  --mapping-suffix MAPPING_SUFFIX
                        mapping files suffix pattern (default: .rmdup.bam)
  --mapping-aligner-options MAPPING_ALIGNER_OPTIONS
                        additional aligner mapping options (default: )

varbin group::
  --varbin-dir VARBIN_DIR
                        varbin working directory (default: None)
  --varbin-suffix VARBIN_SUFFIX
                        varbin files suffix pattern (default: .varbin.txt)
```

### Use of `scclust` subcommand

To list options available for `scclust` subcommand use:

```bash
sgains-tools scclust -h
usage: sgains-tools scclust [-h] [--bins-count BINS_COUNT]
                            [--bins-dir BINS_DIR] [--bins-file BINS_FILE]
                            [--varbin-dir VARBIN_DIR]
                            [--varbin-suffix VARBIN_SUFFIX]
                            [--scclust-case SCCLUST_CASE]
                            [--scclust-dir SCCLUST_DIR]
                            [--scclust-cytoband-file SCCLUST_CYTOBAND_FILE]
                            [--scclust-nsim SCCLUST_NSIM]
                            [--scclust-sharemin SCCLUST_SHAREMIN]
                            [--scclust-fdrthres SCCLUST_FDRTHRES]
                            [--scclust-nshare SCCLUST_NSHARE]
                            [--scclust-climbtoshare SCCLUST_CLIMBTOSHARE]

optional arguments:
  -h, --help            show this help message and exit

bins group::
  --bins-count BINS_COUNT
                        number of bins (default: 10000)
  --bins-dir BINS_DIR   bins working directory (default: None)
  --bins-file BINS_FILE
                        bins boundaries filename (default:
                        bins_boundaries.txt)

varbin group::
  --varbin-dir VARBIN_DIR
                        varbin working directory (default: None)
  --varbin-suffix VARBIN_SUFFIX
                        varbin files suffix pattern (default: .varbin.txt)

scclust group::
  --scclust-case SCCLUST_CASE
                        SCclust case name (default: None)
  --scclust-dir SCCLUST_DIR
                        SCclust working directory (default: None)
  --scclust-cytoband-file SCCLUST_CYTOBAND_FILE
                        location of cyto band description file (default: None)
  --scclust-nsim SCCLUST_NSIM
                        SCclust number of simulations (default: 150)
  --scclust-sharemin SCCLUST_SHAREMIN
                        SCclust sharemin parameter (default: 0.85)
  --scclust-fdrthres SCCLUST_FDRTHRES
                        SCclust fdrthres parameter (default: -3)
  --scclust-nshare SCCLUST_NSHARE
                        SCclust nshare parameter (default: 4)
  --scclust-climbtoshare SCCLUST_CLIMBTOSHARE
                        SCclust climbtoshare parameter (default: 5)
```

## Configure the *s-GAINS* pipeline

An example *s-GAINS* pipeline configuration:

```bash
aligner:
    aligner_name: hisat2

genome:
    genome_version: hg38
    genome_pristine_dir: hg38_pristine
    genome_dir: hg38
    genomeindex_prefix: genomeindex

mappable_regions:
    mappable_read_length: 50
    mappable_dir: hg38_R50
    mappable_file: hisat2_hg38_R50_mappable_regions.txt
    mappable_aligner_options: ""
  
bins:
    bins_count: 20000
    bins_dir: hg38_R50_B20k
    bins_file: hg38_R50_B20k_bins_boundaries.txt

reads:
    reads_dir: navin_T10
    reads_suffix: ".fastq.gz"
    

mapping:
    mapping_dir: navin_T10_hisat2/mapping
    mapping_suffix: ".rmdup.bam"
    mapping_aligner_options: "-3 0 -5 38"

varbin:
    varbin_dir: navin_T10_hisat2/varbin
    varbin_suffix: ".varbin.r50_20k.txt"


scclust:
    scclust_case: "nyu007_hisat2"
    scclust_dir: "navin_T10_hisat2/scclust"
    scclust_cytoband_file: cytoBand-hg38.txt
    scclust_nsim: 150
    scclust_sharemin: 0.85
    scclust_fdrthres: -3
    scclust_nshare: 4
    scclust_climbtoshare: 5
```

Each section of this configuration file corresponds to the relevant `s-GAINS` tool
subcommand and sets values for the options of the subcommand.

The options passed from the command line override the options specified in the
configuration file.

To pass configuration file to `sgains-tools` you should use `-c` or `--config` 
option. For example, if you want to use the config file for `mapping` subcommand
you should use:

```bash
sgains-tools -c sgains-hisat2-navin-T10.yml mapping -h
```

Note that the default values for various parameters of `mapping` subcommand would
be filled from the corresponding values specified into the configuration file:

```bash
usage: sgains-tools mapping [-h] [--aligner-name ALIGNER_NAME]
                            [--genome-version GENOME_VERSION]
                            [--genome-pristine-dir GENOME_PRISTINE_DIR]
                            [--genome-dir GENOME_DIR]
                            [--genomeindex-prefix GENOMEINDEX_PREFIX]
                            [--reads-dir READS_DIR]
                            [--reads-suffix READS_SUFFIX]
                            [--mapping-dir MAPPING_DIR]
                            [--mapping-suffix MAPPING_SUFFIX]
                            [--mapping-aligner-options MAPPING_ALIGNER_OPTIONS]

optional arguments:
  -h, --help            show this help message and exit

aligner group::
  --aligner-name ALIGNER_NAME
                        aligner to use in sGAINS subcommands (default: hisat2)

genome group::
  --genome-version GENOME_VERSION
                        version of reference genome to use (default: hg38)
  --genome-pristine-dir GENOME_PRISTINE_DIR
                        directory where clean copy of reference genome is
                        located (default: /data/lubo/single-
                        cell/test_data/hg38_pristine)
  --genome-dir GENOME_DIR
                        genome index working directory (default:
                        /data/lubo/single-cell/test_data/hg38)
  --genomeindex-prefix GENOMEINDEX_PREFIX
                        genome index prefix (default: genomeindex)

reads group::
  --reads-dir READS_DIR
                        data directory where sequencing reads are located
                        (default: /data/lubo/single-cell/test_data/navin_T10)
  --reads-suffix READS_SUFFIX
                        reads files suffix pattern (default: .fastq.gz)

mapping group::
  --mapping-dir MAPPING_DIR
                        data directory where mapping files are located
                        (default: /data/lubo/single-
                        cell/test_data/navin_T10_hisat2/mapping)
  --mapping-suffix MAPPING_SUFFIX
                        mapping files suffix pattern (default: .rmdup.bam)
  --mapping-aligner-options MAPPING_ALIGNER_OPTIONS
                        additional aligner mapping options (default: -3 0 -5
                        38)

```