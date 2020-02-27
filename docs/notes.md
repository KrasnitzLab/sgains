# Notes

## Hisat2 mapping

```bash
gunzip -c /mnt/wigclust18/data/safe/chorbadj/single_cell_demo/demo_data/cor002/reads/CJA2289.fastq.gz | \
    hisat2 -X 1500 --no-spliced-alignment -3 0 -5 38 -x /mnt/wigclust18/data/safe/chorbadj/single_cell_demo/demo_data/hisat2_hg19/hg19/genomeindex -U - 2> \
    /mnt/wigclust18/data/safe/chorbadj/single_cell_demo/demo_data/hisat2_hg19/cor002/mapping/CJA2289.aligner_report.log | \
    samtools view -bu -q 30 -F 0xff00 | \
    samtools sort | \
    samtools rmdup -s - - | \
    samtools view -b -o /mnt/wigclust18/data/safe/chorbadj/single_cell_demo/demo_data/hisat2_hg19/cor002/mapping/CJA2289.rmdup.bam

samtools index /mnt/wigclust18/data/safe/chorbadj/single_cell_demo/demo_data/hisat2_hg19/cor002/mapping/CJA2289.rmdup.bam
```

## Timing

### Preparation HG38

* mappable regions:

    ```bash
    real    334m34.253s
    user    38m40.777s
    sys     7m44.531s
    ```

* bins:

    ```bash
    real    6m12.246s
    user    0m53.830s
    sys     0m12.830s
    ```

### COR002

* mapping:

    ```bash
    real    108m19.041s
    user    12m10.774s
    sys     2m14.386s
    ```

* varbin:

    ```bash
    real    68m39.065s
    user    8m46.804s
    sys     1m44.197s
    ```

* scclust:

    ```bash
    real    45m32.275s
    user    68m11.955s
    sys     57m24.621s
    ```

### COR003

* mapping

    ```bash
    real    129m21.416s
    user    16m14.520s
    sys     2m59.067s
    ```

* varbin:

    ```bash
    real    77m3.188s
    user    9m51.669s
    sys     1m45.681s
    ```

* scclust:

    ```bash
    real    52m58.345s
    user    80m6.625s
    sys     73m37.284s
    ```
