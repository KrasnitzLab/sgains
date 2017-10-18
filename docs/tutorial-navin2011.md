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

* Download archive file `chromFa.tar.gz` and untar it into a separate directory:

    ```
    mkdir hg19_pristine
    cd hg19_pristine
    wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
    tar zxvf chromFa.tar.gz
    ```
* Go back to the data directory and create a `hg19` subdirectory where the 
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
reference genome. Building this index is computationally intensive process and
could take several hours of CPU time.

### Preparation of uniquely mappable regions

* The next step is to generate uniquely mappable regions, i.e. contiguous regions wherein 
all reads of a given length are unique in the genome. To this end you can use 
`mappalbe-regions` subcommand of `s-GAINS` pipeline. You need to specify the directory, 
where a working copy of reference genome is located and the read length to be used.

* Create a subdirectory in which to save mappable regions file:

    ```
    mkdir R100
    ```

* Here is an example of invoking the `mappalbe-regions` subcommand with reads of
length 100:

    ```
    sgains.py mappable-regions --genome-dir hg19 \
        --mappable-dir R100 --read-length 100
    ````

* This step is computationaly very intensive and could take days in CPU time.
Consider using `--parallel` option of `sgains.py` command to parallelize the
computation if your computer has a suitable number of cores. For example, on a
workstation with 10 cores you could use 8 cores to compute mappable regions:
    ```
    sgains.py -p 8 mappable-regions --genome-dir hg19 \
        --mappable-dir R100 --read-length 100
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

* Create a subdirectory for storing the bin boundaries file:

    ```
    mkdir R100_50k
    ```

* Run `bins` subcommand to calculate bin boundaries.

    ```
    sgains.py bins \
        --mappable-dir R50 \
        --mappable-regions hg19_R100_mappable_regions.txt \
        --genome-dir hg19 \
        --bins-count 50000 \
        --bins-dir R100_50k \
        --bins-boundaries hg19_R50_B50k_bins_boundaries.txt
    ```
    
    ```
    sgains.py bins \
        --mappable-dir R100 \
        --mappable-regions hg19_R100_mappable_regions.txt \
        --genome-dir hg19 \
        --bins-count 10000 \
        --bins-dir R100_B10k \
        --bins-boundaries hg19_R100_B10k_bins_boundaries.txt
    ```

* To run the command you need to specify:
    * the number of bins you want to calculate
    * a directory for storing the bin boundary file
    * a directory and file name where mappble regions file name is located
    * a directory where a working copy of HG19 is located


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

* Create a subdirectory `Navin2011_T10` for storing the downloaded read file:

    ```
    mkdir Navin2011_T10
    cd Navin2011_T10
    ```

* To download read files for T10 Ductal Carcinoma you can use `fastq-tool`. If
you need a read file for a single sample, you can use:

    ```
    fastq-dump --gzip SRR089402
    ```
This command will download a read file in `fastq` format for a sample with 
accession number *SRR089402*. 

* If you want to download read files for all samples from accession list 
`SRR_Acc_List.txt`, you can use:

    ```
    cat SRR_Acc_List.txt | xargs fastq-dump --gzip
    ```

---

**Please note that the last command will download about 50Gb of data and will store
about 100Gb of data on disk (cache and actual reads).**

---


### Process downloaded data

* To process downloaded data you can use the `process` subcommand:

    ```
    sgains.py -p 8 process --genome-dir hg19 \
        --bins-dir R50_50k \
        --mapping-bowtie-opts "-S -t -m 1 --best --strata" \
        --reads-dir Navin2011_T10 \
        --study-name "navin2011_t10" -o Navin2011_T10_Results
    ```

* Note that your configuration file contains values for most of the 
pipeline parameters. So, once you have a configuration file you can skip most 
of the parameters and use:

    ```
    sgains.py -p 8 process -o Navin2011_T10_Results
    ```
