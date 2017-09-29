
## Download data for T10 Ductal Carcinoma

Go to 
```
https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SAMN00014736
```

and use *Accession List* button to download file `SRR_Acc_List.txt` containing
all samples accession numbers for this experiment.

To download the sample reads you need to use *SRA Toolkit*. SRA Toolkit is
available in Anaconda `bioconda` channel. You can install it using your Anaconda 
environment `sgains`:

```
conda install sra-tools
```

To downalod samples reads for T10 Ductal Carcinoma you can use `fastq-tool`. If
you need single sample read you can use:

```
fastq-dump --gzip SRR089402
```

This command will download sample reads in `fastq` format for sample with 
accession number *SRR089402*. 

If you want to download all samples from accession list `SRR_Acc_List.txt` you
can use:

```
cat SRR_Acc_List.txt | xargs fastq-dump --gzip
```

Please note, that this command will download about 50Gb of data and will store
about 100Gb of data on disk (cache and actual reads).

