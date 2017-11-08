# Installation of `s-GAINS` development environment on Mac OS X

## Install R

Download *R* installation pacakges from 
[https://cran.r-project.org/bin/macosx/](https://cran.r-project.org/bin/macosx/)
and follow instructions for installing *R* on your Mac.

## Install Anaconda 3.6

Download *Anaconda* installation package for *Python 3.6* from 
[https://www.anaconda.com/download/#macos](https://www.anaconda.com/download/#macos)
and follow instructions for installtion *Anaconda* on your Mac.

## Install required python packages

Create a `conda` environment for working with *s-GAINS* pipeline:

```
conda create -n sgains3
```

Activate this environment:
```
source activate sgains3
```

and install following packages into this working environment:

```
conda install pandas
conda config --add channels bioconda
conda install samtools bcftools biopython pysam

conda install -c conda-forge perl=5.22.0
conda install bowtie=1.2.1.1

pip install python-box termcolor PyYAML pytest pytest-asyncio
```

## Get the `s-GAINS` source code 

You can use git to clone the source code of `sgains` pipeline:

```
git clone https://github.com/KrasnitzLab/sgains.git
```

Alternatively you can download the source code from
[https://github.com/KrasnitzLab/sgains](https://github.com/KrasnitzLab/sgains)
and unarhive the downloaded ZIP file.

Let assume that the source code is located into a subdirectory `sgains`. 
Directory structure is as follows:

```
.
└── sgains
    ├── LICENSE
    ├── README.md
    ├── data
    ├── docs
    ├── scpipe
    ├── scripts
    ├── setenv.sh
    ├── sgains.yml
    └── tools
```

## Install additional R packages

The `s-GAINS` pipeline depends on `DNAcopy` and `TBEST`. Additionally you should
install `SCclust` clust package.

Directory `sgains/scripts` contains and `R` script, that installs all required
*R* packages. Enter into this directory and execute following command:

```
cd sgains/scripts
Rscript setup.R
```

## Update your shell environment

To work with `s-GAINS` command line tool `sgains.py` you need to make some changes
in your shell environment. 

You need to add `sgains/tools` directory
```
sgains/
└── tools
```
to the `PATH` variable and also you need to add `sgains/scpipe` subdirectory
```
sgains/
└── scpipe
```

to the `PYTHONPATH` variable.

To make this changes permanent you can edit `$HOME/.bash_profile`. Let assume
that the source code of `s-GAINS` pipeline is located inside your home directory
`$HOME/sgains`. Then you can add following lines to your `.bash_profile` file:

```
...
export PATH=$HOME/sgains/tools:$PATH
export PYTHONPATH=$HOME/sgains/scpipe:$PYTHONPATH
```

## Check your environment

After these steps you should be able to call `sgains.py` command from your
terminal:

    ```
    sgains.py --help
    usage: sgains.py [-h] [-v] [-c path] [-n] [--force] [--parallel PARALLEL]
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
      --parallel PARALLEL, -p PARALLEL
                            number of task to run in parallel
    
    subcommands:
      {process,prepare,genomeindex,mappable-regions,bins,mapping,varbin,segment}
        process             combines mapping, varbin and segment subcommands into
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
        segment             segments bin counts and prepares the SCGV input data
    ```
