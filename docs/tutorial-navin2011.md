## `s-GAINS` pipeline prepration

Create a directory where to place all data files you plan to process with 
`s-GAINS` pipeline:

```
mkdir data
cd data
```

All instructions bellow asume that you are working inside your data directory.

### Preparation of genome index

* First you need a copy of human reference genome version `hg19`. To download
it you can go to [UCSC Genome Browser](https://genome.ucsc.edu/), locate the 
downloads section and find full data set for *GRCh37/hg19* version of human 
reference genome. 

* Download archive file `chromFa.tar.gz` and untar it into separate directory:

    ```
    mkdir hg19_pristine
    cd hg19_pristine
    wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
    tar zxvf chromFa.tar.gz
    ```
* Go back to data directory and create a `hg19` subdirectory where the 
pipeline will place a modified version of `hg19` reference genome:

    ```
    mkdir hg19
    ```

* Run `genomeindex` subcommand to copy and modify `hg19` reference genome from
your `hg19_pristine` copy into working `hg19` subdirectory:

    ```
    sgains.py genomeindex --genome-pristine hg19_pristine --genome-dir hg19
    ```

* This step will use `bowtie-build` command produce bowtie index of the *hg19* 
reference genome. Building this index is computationally intesive process and
could take several hours of CPU time.

### Preparation of uniquely mappable regions

* The next step is to generate uniquely mappable regions with given lenght. To 
this end you can use `mappalbe-regions` subcommand of `s-GAINS` pipeline. You
need to specify the directory, where working copy of reference genome is located
and the length of mappalbe regions you need.

* Create subdirectory where to save mappable regions file:

    ```
    mkdir R50
    ```

* Example invocation of `mappalbe-regions` subcommand specifying mappable regions
with length 100 base pairs is:

    ```
    sgains.py mappable-regions --genome-dir hg19 \
        --mappable-dir R50 --read-length 50
    ````

* This step is computationaly very intensive and could take days in CPU time.
Consider using `--parallel` option of `sgains.py` command to parallelize the
computations if your computer have suitable number of cores. For example, on a
workstation with 10 cores you could for example use 8 cores for processing 
mappable regions:
    ```
    sgains.py -p 8 mappable-regions --genome-dir hg19 \
        --mappable-dir R50 --read-length 50
    ```

* Alternatively you can download precomputed mappable regions file from `s-GAINS`
pipeline releases at [https://github.com/KrasnitzLab/sgains/releases](https://github.com/KrasnitzLab/sgains/releases). 
For example you can download mappable regions with length 100 base pairs for HG19
reference genome from [https://github.com/KrasnitzLab/sgains/releases/download/1.0_beta1/hg19_R100_mappable_regions.txt.gz].

    ```
    mkdir R100d
    cd R100d/
    wget -c https://github.com/KrasnitzLab/sgains/releases/download/1.0_beta1/hg19_R100_mappable_regions.txt.gz
    gunzip hg19_R100_mappable_regions.txt.gz
    ```

### Calculation of bins boundaries

* Create a subdirectory where to store bins boundaries file:

    ```
    mkdir R50_50k
    ```

* Run `bins` subcommand to calculate bins boundaries.

    ```
    sgains.py bins --mappable-dir R50 \
        --genome-dir hg19 \
        --bins-dir R50_50k --bins-count 50000 \
        --bins-boundaries hg19_R50_50k_bins_boundaries.txt
    ```

* To run the command you need to specify:
    * the number of bins you want to calculate
    * directory where to store the bins file
    * directory and file name where mappble regions file name is located
    * directory where working copy of HG19 is located


## Processing data with `s-GAINS` pipeline

To demonstrate the usage of `s-GAINS` pipeline we are going to use data from
*Tumour evolution inferred by single-cell sequencing* paper ([https://www.ncbi.nlm.nih.gov/pubmed/21399628](https://www.ncbi.nlm.nih.gov/pubmed/21399628))

### Configure the pipeline

Since the pipeline has many parameters you can create a configuration file, that
sets values for most of parameters used by the pipeline.

The configuration file is in YAML format and has the following structure:

```
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
    bins_count: 50000
    bins_dir: R50_50k
    bins_boundaries: hg19_R50_50k_bins_boundaries.txt

mapping:
    reads_dir: Navin2011_T10
    reads_suffix: ".fastq.gz"
    mapping_dir: mappings
    mapping_suffix: ".rmdup.bam"
    mapping_bowtie_opts: "-S -t -m 1 --best --strata"

varbin:
    varbin_dir: varbin
    varbin_suffix: ".varbin.50k.txt"

segment:
    segment_dir: segment
    study_name: "navin2011_T10"

```

### Download data for T10 Ductal Carcinoma

* Go to [https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SAMN00014736](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SAMN00014736) 
and use *Accession List* button to download file `SRR_Acc_List.txt` containing
all samples accession numbers for this experiment.

* To download the sample reads you need to use *SRA Toolkit*. SRA Toolkit is
available in Anaconda `bioconda` channel. You can install it using your Anaconda 
environment `sgains`:

    ```
    conda install sra-tools
    ```

* Create a subdirectory `Navin2011_T10` where to store downloaded samples file:

    ```
    mkdir Navin2011_T10
    cd Navin2011_T10
    ```

* To downalod samples reads for T10 Ductal Carcinoma you can use `fastq-tool`. If
you need single sample read you can use:

    ```
    fastq-dump --gzip SRR089402
    ```
This command will download sample reads in `fastq` format for sample with 
accession number *SRR089402*. 

* If you want to download all samples from accession list `SRR_Acc_List.txt` you
can use:

    ```
    cat SRR_Acc_List.txt | xargs fastq-dump --gzip
    ```

---

**Please note, that last command will download about 50Gb of data and will store
about 100Gb of data on disk (cache and actual reads).**

---


### Process downloaded data

* To process downloaded data you can use `process` subcommand:

    ```
    sgains.py -p 8 process --genome-dir hg19 \
        --bins-dir R50_50k \
        --mapping-bowtie-opts "-S -t -m 1 --best --strata" \
        --reads-dir Navin2011_T10 \
        --study-name "navin2011_t10" -o Navin2011_T10_Results
    ```

* Note that your configuration file contains values for most of the 
pipeline parameters. So once you have a configuration file you can skip most 
of the parameters and use:

    ```
    sgains.py -p 8 process -o Navin2011_T10_Results
    ```
