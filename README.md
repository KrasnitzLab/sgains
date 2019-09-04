# Sparse Genomic Analysis of Individual Nuclei by Sequencing (s-GAINS)

[![DOI](https://zenodo.org/badge/98500084.svg)](https://zenodo.org/badge/latestdoi/98500084)

This document describes how to setup `s-GAINS` pipeline tool and its basic command.

Short tutorial on how to use this tool could be found in
[Example usage of `sGAINS` pipeline](docs/tutorial-navin2011.md)

## Anaconda environment setup

### Install Anaconda

* Go to anaconda web site
    [https://www.continuum.io/downloads](https://www.continuum.io/downloads)
    and download the latest anaconda installer for your operating system.

* *s-GAINS* supports *Python 3.6* or greater so you need to choose an
    appropriate installer. Note also that since *s-GAINS* uses *bioconda*
    channel the supported operating systems are only those supported for
    *bioconda* (at the time of this writing these are Linux and Mac OS X).

* Install anaconda into suitable place on your local machine following
    instructions from
    [https://docs.continuum.io/anaconda/install](https://docs.continuum.io/anaconda/install)

### Create `sgains` Anaconda environment

* After installing and activating *Anaconda* you need to create an environment to
    use with `sgains` pipeline. To this end you need to use:

    ```bash
    conda create -n sgains3
    source activate sgains3
    ```

### Install `sgains` anaconda package

* *sGAINS* tools are distributed as a conda package through `krasnitzlab`
    Annaconda channel. So to install *sGAINS* tools use:

    ```bash
    conda install -c krasnitzlab -c bioconda -c conda-forge sgains
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
                 {process,prepare,genomeindex,mappable-regions,bins,mapping,varbin,scclust}
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
  --parallel PARALLEL, -p PARALLEL
                        number of task to run in parallel

subcommands:
  {process,prepare,genomeindex,mappable-regions,bins,mapping,varbin,scclust}
    process             combines mapping, varbin and scclust subcommands into
                        single command
    prepare             combines all preparation steps (genomeindex,
                        mappable_regions, bins) into single command
    genomeindex         builds appropriate bowtie index for the reference
                        genome
    mappable-regions    finds all mappable regions in specified genome
    bins                calculates all bins boundaries for specified bins
                        count and read length
    mapping             performs mapping of cell reads to reference genome
    varbin              applies varbin algorithm to count read mappings in
                        each bin
    scclust             segmentation and clustering based bin counts and
                        preparation of the SCGV input data```
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

## Pipeline preparation

### Usage of `genomeindex` subcommand

The `genomeindex` subcommand builds the bowtie index for the reference genome. To
list the available options use:

```bash
sgains-tools genomeindex -h
usage: sgains-tools genomeindex [-h] [--genome-index GENOME_INDEX]
                             [--genome-dir GENOME_DIR]
                             [--genome-version GENOME_VERSION]
                             [--genome-pristine GENOME_PRISTINE]

optional arguments:
  -h, --help            show this help message and exit

genome index options:
  --genome-index GENOME_INDEX, -G GENOME_INDEX
                        genome index name (default: genomeindex)
  --genome-dir GENOME_DIR
                        genome index directory (default: hg19)
  --genome-version GENOME_VERSION
                        version of reference genome in use (supports only
                        hg19) (default: hg19)
  --genome-pristine GENOME_PRISTINE
                        directory where clean copy of reference genome is
                        located (default: hg19_pristine)
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
usage: sgains-tools mappable-regions [-h] [--mappable-dir MAPPABLE_DIR]
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
                        (default: R50)
  --mappable-regions MAPPABLE_REGIONS, -M MAPPABLE_REGIONS
                        filename where mappable regions are stored (default:
                        hg19_R100_mappable_regions.txt)
  --read-length LENGTH, -l LENGTH
                        read length to use for generation of mappable regions
                        (default: 50)
  --bowtie-opts BOWTIE_OPTS
                        additional bowtie options (default: )

genome index options:
  --genome-index GENOME_INDEX, -G GENOME_INDEX
                        genome index name (default: genomeindex)
  --genome-dir GENOME_DIR
                        genome index directory (default: hg19)
  --genome-version GENOME_VERSION
                        version of reference genome in use (supports only
                        hg19) (default: hg19)
```

### Usage of `bins` subcommand

The `bins` subcommand computes the bins boudaries.

To list options available for `bins` subcommand use:

```bash
sgains-tools bins -h
usage: sgains-tools bins [-h] [--mappable-dir MAPPABLE_DIR]
                      [--mappable-regions MAPPABLE_REGIONS]
                      [--bins-boundaries BINS_BOUNDARIES]
                      [--bins-dir BINS_DIR] [--bins-count BINS_COUNT]
                      [--genome-index GENOME_INDEX] [--genome-dir GENOME_DIR]
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

bins boundaries:
  --bins-boundaries BINS_BOUNDARIES, -B BINS_BOUNDARIES
                        bins boundaries filename (default:
                        bins_boundaries.tst)
  --bins-dir BINS_DIR   bins working directory (default: data/R100_B10k)
  --bins-count BINS_COUNT, -C BINS_COUNT
                        number of bins (default: 10000)

genome index options:
  --genome-index GENOME_INDEX, -G GENOME_INDEX
                        genome index name (default: genomeindex)
  --genome-dir GENOME_DIR
                        genome index directory (default: data/hg19)
  --genome-version GENOME_VERSION
                        version of reference genome in use (supports only
                        hg19) (default: hg19)
```

## Processing sequence data

### Use of `process` subcommand

---

**Please note, that to use `process` subcommands
(`mapping`, `varbin`, `scclust` and `process`)
you need go through preparation steps.**

---

To list the options available for `process` subcommand use:

```bash
sgains-tools process -h
usage: sgains-tools process [-h] [--reads-dir READS_DIR]
                         [--reads-suffix READS_SUFFIX]
                         [--mapping-bowtie-opts MAPPING_BOWTIE_OPTS]
                         [--output-dir OUTPUT_DIR] [--case-name CASE_NAME]
                         [--genome-index GENOME_INDEX]
                         [--genome-dir GENOME_DIR]
                         [--genome-version GENOME_VERSION]
                         [--bins-boundaries BINS_BOUNDARIES]
                         [--bins-dir BINS_DIR]

optional arguments:
  -h, --help            show this help message and exit
  --mapping-bowtie-opts MAPPING_BOWTIE_OPTS
                        bowtie mapping options (default: -S -t -m 1 --best
                        --strata --chunkmbs 256)

sequencing reads options:
  --reads-dir READS_DIR, -R READS_DIR
                        data directory where sequencing reads are located
                        (default: SRA)
  --reads-suffix READS_SUFFIX
                        reads files suffix pattern (default: .fastq.gz)

process output options:
  --output-dir OUTPUT_DIR, -o OUTPUT_DIR
                        output directory (default: scgv)
  --case-name CASE_NAME
                        case name (default: navin_T10)

genome index options:
  --genome-index GENOME_INDEX, -G GENOME_INDEX
                        genome index name (default: genomeindex)
  --genome-dir GENOME_DIR
                        genome index directory (default: hg19)
  --genome-version GENOME_VERSION
                        version of reference genome in use (supports only
                        hg19) (default: hg19)

bins boundaries:
  --bins-boundaries BINS_BOUNDARIES, -B BINS_BOUNDARIES
                        bins boundaries filename (default:
                        hg19_R50_B20k_bins_boundaries.txt)
  --bins-dir BINS_DIR   bins working directory (default: R50_B20k)
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
usage: sgains-tools mapping [-h] [--reads-dir READS_DIR]
                         [--reads-suffix READS_SUFFIX]
                         [--mapping-dir MAPPING_DIR]
                         [--mapping-suffix MAPPING_SUFFIX]
                         [--mapping-bowtie-opts MAPPING_BOWTIE_OPTS]
                         [--genome-index GENOME_INDEX]
                         [--genome-dir GENOME_DIR]
                         [--genome-version GENOME_VERSION]

optional arguments:
  -h, --help            show this help message and exit
  --mapping-bowtie-opts MAPPING_BOWTIE_OPTS
                        bowtie mapping options (default: -S -t -m 1 --best
                        --strata)

sequencing reads options:
  --reads-dir READS_DIR, -R READS_DIR
                        data directory where sequencing reads are located
                        (default: reads)
  --reads-suffix READS_SUFFIX
                        reads files suffix pattern (default: .fastq.gz)

mapping files options:
  --mapping-dir MAPPING_DIR, -M MAPPING_DIR
                        data directory where mapping files are located
                        (default: mappings)
  --mapping-suffix MAPPING_SUFFIX
                        mapping files suffix pattern (default: .rmdup.bam)

genome index options:
  --genome-index GENOME_INDEX, -G GENOME_INDEX
                        genome index name (default: genomeindex)
  --genome-dir GENOME_DIR
                        genome index directory (default: ../../hg19)
  --genome-version GENOME_VERSION
                        version of reference genome in use (supports only
                        hg19) (default: hg19)
```

### Use of `varbin` subcommand

To list the options available for `varbin` subcommand use:

```bash
sgains-tools varbin -h
usage: sgains-tools varbin [-h] [--mapping-dir MAPPING_DIR]
                        [--mapping-suffix MAPPING_SUFFIX]
                        [--varbin-dir VARBIN_DIR]
                        [--varbin-suffix VARBIN_SUFFIX]
                        [--bins-boundaries BINS_BOUNDARIES]
                        [--bins-dir BINS_DIR]

optional arguments:
  -h, --help            show this help message and exit

mapping files options:
  --mapping-dir MAPPING_DIR, -M MAPPING_DIR
                        data directory where mapping files are located
                        (default: mappings)
  --mapping-suffix MAPPING_SUFFIX
                        mapping files suffix pattern (default: .rmdup.bam)

varbin options:
  --varbin-dir VARBIN_DIR, -V VARBIN_DIR
                        varbin directory (default: varbin)
  --varbin-suffix VARBIN_SUFFIX
                        varbin files suffix pattern (default: .varbin.10k.txt)

bins boundaries:
  --bins-boundaries BINS_BOUNDARIES, -B BINS_BOUNDARIES
                        bins boundaries filename (default:
                        hg19_R50_B50k_bins_boundaries.txt)
  --bins-dir BINS_DIR   bins working directory (default: ../../R50_B50k)
```

### Use of `scclust` subcommand

To list options available for `scclust` subcommand use:

```bash
sgains-tools scclust -h
usage: sgains-tools scclust [-h] [--varbin-dir VARBIN_DIR]
                         [--varbin-suffix VARBIN_SUFFIX]
                         [--scclust-dir SCCLUST_DIR] [--case-name CASE_NAME]
                         [--bins-boundaries BINS_BOUNDARIES]
                         [--bins-dir BINS_DIR]

optional arguments:
  -h, --help            show this help message and exit

varbin options:
  --varbin-dir VARBIN_DIR, -V VARBIN_DIR
                        varbin directory (default: varbin)
  --varbin-suffix VARBIN_SUFFIX
                        varbin files suffix pattern (default: .varbin.20k.txt)

SCclust options:
  --scclust-dir SCCLUST_DIR, -S SCCLUST_DIR
                        SCGV directory (default: scclust)
  --case-name CASE_NAME
                        case name (default: navin_T10)

bins boundaries:
  --bins-boundaries BINS_BOUNDARIES, -B BINS_BOUNDARIES
                        bins boundaries filename (default:
                        hg19_R50_B20k_bins_boundaries.txt)
  --bins-dir BINS_DIR   bins working directory (default: R50_B20k)
```

## Configure the *s-GAINS* pipeline

An example *s-GAINS* pipeline configuration:

```bash
genome:
    version: hg19
    data_dir: hg19_pristine
    work_dir: hg19
    genome_index: genomeindex

mappable_regions:
    length: 50
    work_dir: R50
    bowtie_opts: ""
  
bins:
    bins_count: 20000
    bins_dir: R50_B20k
    bins_boundaries: hg19_R50_B20k_bins_boundaries.txt

mapping:
    reads_dir: SRA
    reads_suffix: ".fastq.gz"
    mapping_dir: mapping
    mapping_suffix: ".rmdup.bam"
    mapping_bowtie_opts: "-S -t -m 1 --best --strata --chunkmbs 256"

varbin:
    varbin_dir: varbin
    varbin_suffix: ".varbin.20k.txt"

scclust:
    case_name: "navin_T10"
    scgv_dir: scgv
    cytoband: hg19/cytoBand.txt
    nsim: 150
    sharemin: 0.85
```

Each section of this configuration file corresponds to the relevant `s-GAINS` tool
subcommand and sets values for the options of the subcommand.

The options passed from the command line override the options specified in the
configuration file.
