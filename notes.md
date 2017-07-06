

## Install tools

```
conda config --add channels bioconda
conda install samtools bcftools
conda install wgsim
conda install trimmomatic
conda install biopython pysam

conda install fastqc

```

## Nature protocols paper

```
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5069701/
```

## Other papers

```
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4504184/
```


## Sam tools tutorial

```
http://biobits.org/samtools_primer.html
```
Interesting to check this:
```
https://wikis.utexas.edu/display/bioiteam/Bioinformatics+Courses+and+Content
https://wikis.utexas.edu/display/bioiteam/Genome+Variant+Analysis+Course+2017
https://wikis.utexas.edu/display/bioiteam/Introduction+to+RNA+Seq+Course+2017
```

More tutorials:
```
https://github.com/davetang/learning_bam_file

```


## Build the HG19 bowtie index

## Generate reads

```
cd data/hg19
python generate.reads.py 100 chr1 > ../readsim/chr01.reads.fa 
```

The script `generate.reads.py` uses about 10G RAM. 


## Align reads with `bowtie`

```
cd data/hg19
bowtie -S -t -v 0 -m 1 -f genomeindex ../readsim/chr01.reads.fa > ../readsim/chr01.reads.mapping.sam
```

## Mappable regions

```
python mappable.regions.py < data/readsim/chr01.reads.mapping.sam > data/readsim/chr01.stdout.txt 2> data/readsim/chr01.stderr.txt
```

## All mappable regions

```

```

## Sort mappable regions

```
sort -k 1,1 -k 2,2n data/readsim/chr01.stdout.txt > mappable.regions.sorted.txt
```


## Experiment with my tool chain

```
generate_mappable_regions.py -c scpipe.yml -l 100 -C chr1 > data/chr1.stdout.txt
```

```
generate_reads.py -c /mnt1/lubo/Work/single_cell_pipeline/scpipe.yml -C chr1 -l 100 | \
    bowtie2 -t -x genomeindex -f - | \
    mappable_regions.py -c /mnt1/lubo/Work/single_cell_pipeline/scpipe.yml
```

The original tool chain whould be something like this:

```
generate_reads_jc.py 100 hg19 | bowtie2 -t -x genomeindex -f - | mappable_regions_jc.py > chr1.orig.stdout.txt
```

## Time consumed to generate mappable regions

```
date && generate_mappable_regions.py -l 50 -c scpipe.yml -o data/readsim/all.50.mappable.txt && date
Fri Jun 16 14:14:26 EEST 2017
...
...
Sat Jun 17 04:20:49 EEST 2017
```

## Difference in results when writing reads

When we write reads with *biopython* `SeqIO` module  and when we write reads
by hand, the results of mappings is different. Not sure how bowtie works.



## FASTQ quality scores

```
https://en.wikipedia.org/wiki/FASTQ_format
```


## Illumina adapter trimming

* Trimmomatic: a flexible trimmer for Illumina sequence data
```
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4103590/
```

* Trimmomatic manual:
```
http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf
```

* Question: How To Best Deal With Adapter Contamination (Illumina)?
```
https://www.biostars.org/p/3461/
```

* Adapter and quality trimming of illumina data
```
http://www.ark-genomics.org/events-online-training-eu-training-course/adapter-and-quality-trimming-illumina-data
```

## Generation of mappable regions - 100bp and 50bp

```
date && generate_mappable_regions.py -c scpipe10k100.yml -l 100 && date
Fri Jun 23 14:56:20 EEST 2017
Time loading forward index: 00:00:01
Time for 0-mismatch search: 32:57:10
# reads processed: 3095675036
# reads with at least one reported alignment: 2761401626 (89.20%)
# reads that failed to align: 234379429 (7.57%)
# reads with alignments suppressed due to -m: 99893981 (3.23%)
Reported 2761401626 alignments to 1 output stream(s)
Time searching: 32:57:11
Overall time: 32:57:11
Sat Jun 24 23:53:31 EEST 2017
```

```
date && generate_mappable_regions.py -c scpipe10k50.yml -l 50 && date                                                                                       
Mon Jun 26 16:55:47 EEST 2017
Time loading forward index: 00:00:01
Time for 0-mismatch search: 32:23:16
# reads processed: 3095676236
# reads with at least one reported alignment: 2625057481 (84.80%)
# reads that failed to align: 234364827 (7.57%)
# reads with alignments suppressed due to -m: 236253928 (7.63%)
Reported 2625057481 alignments to 1 output stream(s)
Time searching: 32:23:17
Overall time: 32:23:17
Wed Jun 28 01:19:12 EEST 2017
```


## Timing of mappings generation

```
date && generate_mappings.py -c scpipe_tests.yml -l 100 -C chr4 | samtools view -o test_mappings_chr4.bam - && date                                         
Tue Jun 27 11:16:41 EEST 2017
Time loading forward index: 00:00:28
Time for 0-mismatch search: 02:33:22
# reads processed: 191154177
# reads with at least one reported alignment: 183818120 (96.16%)
# reads that failed to align: 3493590 (1.83%)
# reads with alignments suppressed due to -m: 3842467 (2.01%)
Reported 183818120 alignments to 1 output stream(s)
Time searching: 02:33:57
Overall time: 02:33:57
Tue Jun 27 13:50:44 EEST 2017
```

```
date && generate_mappings.py -c scpipe_tests.yml -l 100 -C chrY | samtools view -o test_mappings_chrY.bam - && date                                         
Tue Jun 27 11:15:19 EEST 2017
Time loading forward index: 00:00:10
Time for 0-mismatch search: 00:47:40
# reads processed: 59373467
# reads with at least one reported alignment: 17264921 (29.08%)
# reads that failed to align: 33721683 (56.80%)
# reads with alignments suppressed due to -m: 8386863 (14.13%)
Reported 17264921 alignments to 1 output stream(s)
Time searching: 00:47:50
Overall time: 00:47:51
Tue Jun 27 12:03:17 EEST 2017
```

```
date && generate_mappings.py -c scpipe_tests.yml -l 100 -C chrX | samtools view -o test_mappings_chrX.bam - && date                                         
Tue Jun 27 11:14:20 EEST 2017
Time loading forward index: 00:00:02
Time for 0-mismatch search: 02:25:42
# reads processed: 155270461
# reads with at least one reported alignment: 143846327 (92.64%)
# reads that failed to align: 4172079 (2.69%)
# reads with alignments suppressed due to -m: 7252055 (4.67%)
Reported 143846327 alignments to 1 output stream(s)
Time searching: 02:25:45
Overall time: 02:25:45
Tue Jun 27 13:40:13 EEST 2017
```

## Jude's Reads alignment

```
bash commands.sh 
Time loading forward index: 00:00:09
Time loading mirror index: 00:00:09
Seeded quality full-index search: 00:36:16
# reads processed: 5663662
# reads with at least one reported alignment: 1347631 (23.79%)
# reads that failed to align: 3910997 (69.05%)
# reads with alignments suppressed due to -m: 405034 (7.15%)
Reported 1347631 alignments to 1 output stream(s)
Time searching: 00:36:34
Overall time: 00:36:34
Time loading forward index: 00:00:09
Time loading mirror index: 00:00:08
Seeded quality full-index search: 00:19:20
# reads processed: 3805696
# reads with at least one reported alignment: 1385167 (36.40%)
# reads that failed to align: 1947063 (51.16%)
# reads with alignments suppressed due to -m: 473466 (12.44%)
Reported 1385167 alignments to 1 output stream(s)
Time searching: 00:19:37
Overall time: 00:19:37
Time loading forward index: 00:00:04
Time loading mirror index: 00:00:08
Seeded quality full-index search: 00:15:32
# reads processed: 3103879
# reads with at least one reported alignment: 1028255 (33.13%)
# reads that failed to align: 1749477 (56.36%)
# reads with alignments suppressed due to -m: 326147 (10.51%)
Reported 1028255 alignments to 1 output stream(s)
Time searching: 00:15:44
Overall time: 00:15:44
Time loading forward index: 00:00:04
Time loading mirror index: 00:00:06
Warning: Exhausted best-first chunk memory for read TUPAC:142:FC64VJLAAXX:2:108:13313:13867 1:N:0: (patid 2935437); skipping read
Seeded quality full-index search: 00:16:52
# reads processed: 3298881
# reads with at least one reported alignment: 1051553 (31.88%)
# reads that failed to align: 1848631 (56.04%)
# reads with alignments suppressed due to -m: 398697 (12.09%)
Reported 1051553 alignments to 1 output stream(s)
Time searching: 00:17:02
Overall time: 00:17:02
```

Using the other adapter:

```
bash commands2.sh 
Time loading forward index: 00:00:02
Time loading mirror index: 00:00:09
Warning: Exhausted best-first chunk memory for read TUPAC:114:FC64FLDAAXX:3:64:9169:11935 1:N:0: (patid 2427728); skipping read
Seeded quality full-index search: 00:25:38
# reads processed: 4480186
# reads with at least one reported alignment: 1278512 (28.54%)
# reads that failed to align: 2804646 (62.60%)
# reads with alignments suppressed due to -m: 397028 (8.86%)
Reported 1278512 alignments to 1 output stream(s)
Time searching: 00:25:49
Overall time: 00:25:49
Time loading forward index: 00:00:05
Time loading mirror index: 00:00:07
Seeded quality full-index search: 00:23:09
# reads processed: 4697846
# reads with at least one reported alignment: 1785591 (38.01%)
# reads that failed to align: 2406588 (51.23%)
# reads with alignments suppressed due to -m: 505667 (10.76%)
Reported 1785591 alignments to 1 output stream(s)
Time searching: 00:23:21
Overall time: 00:23:21
Time loading forward index: 00:00:04
Time loading mirror index: 00:00:08
Seeded quality full-index search: 00:15:34
# reads processed: 3193189
# reads with at least one reported alignment: 1253946 (39.27%)
# reads that failed to align: 1548138 (48.48%)
# reads with alignments suppressed due to -m: 391105 (12.25%)
Reported 1253946 alignments to 1 output stream(s)
Time searching: 00:15:46
Overall time: 00:15:46
Time loading forward index: 00:00:03
Time loading mirror index: 00:00:07
Seeded quality full-index search: 00:14:15
# reads processed: 2787738
# reads with at least one reported alignment: 1024765 (36.76%)
# reads that failed to align: 1434114 (51.44%)
# reads with alignments suppressed due to -m: 328859 (11.80%)
Reported 1024765 alignments to 1 output stream(s)
Time searching: 00:14:25
Overall time: 00:14:25
Time loading forward index: 00:00:03
Time loading mirror index: 00:00:05
Seeded quality full-index search: 00:20:21
# reads processed: 3441379
# reads with at least one reported alignment: 832358 (24.19%)
# reads that failed to align: 2302985 (66.92%)
# reads with alignments suppressed due to -m: 306036 (8.89%)
Reported 832358 alignments to 1 output stream(s)
Time searching: 00:20:29
Overall time: 00:20:29
Time loading forward index: 00:00:01
Time loading mirror index: 00:00:06
Seeded quality full-index search: 00:36:01
# reads processed: 5673613
# reads with at least one reported alignment: 1350221 (23.80%)
# reads that failed to align: 3914178 (68.99%)
# reads with alignments suppressed due to -m: 409214 (7.21%)
Reported 1350221 alignments to 1 output stream(s)
Time searching: 00:36:08
Overall time: 00:36:08
Time loading forward index: 00:00:06
Time loading mirror index: 00:00:06
Seeded quality full-index search: 00:19:06
# reads processed: 3818801
# reads with at least one reported alignment: 1389024 (36.37%)
# reads that failed to align: 1952788 (51.14%)
# reads with alignments suppressed due to -m: 476989 (12.49%)
Reported 1389024 alignments to 1 output stream(s)
Time searching: 00:19:18
Overall time: 00:19:18
Time loading forward index: 00:00:03
Time loading mirror index: 00:00:07
Seeded quality full-index search: 00:15:27
# reads processed: 3123316
# reads with at least one reported alignment: 1035829 (33.16%)
# reads that failed to align: 1755310 (56.20%)
# reads with alignments suppressed due to -m: 332177 (10.64%)
Reported 1035829 alignments to 1 output stream(s)
Time searching: 00:15:37
Overall time: 00:15:37
Time loading forward index: 00:00:04
Time loading mirror index: 00:00:04
Warning: Exhausted best-first chunk memory for read TUPAC:142:FC64VJLAAXX:2:108:13313:13867 1:N:0: (patid 2941767); skipping read
Seeded quality full-index search: 00:16:44
# reads processed: 3306021
# reads with at least one reported alignment: 1053262 (31.86%)
# reads that failed to align: 1852040 (56.02%)
# reads with alignments suppressed due to -m: 400719 (12.12%)
Reported 1053262 alignments to 1 output stream(s)
Time searching: 00:16:52
Overall time: 00:16:52
```
