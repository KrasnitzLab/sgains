#!/usr/bin/env python

import sys


def main():
    read_length = int(sys.argv[1])
    ref_genome = sys.argv[2]
    if ref_genome == "hg19":
        list_dir = [
            "chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
            "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13",
            "chr14", "chr15", "chr16", "chr17", "chr18",
            "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
    elif ref_genome == "hgdm":
        list_dir = [
            "chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
            "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13",
            "chr14", "chr15", "chr16", "chr17", "chr18",
            "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "dm_chrX",
            "dm_chr2L", "dm_chr2R", "dm_chr3L", "dm_chr3R", "dm_chr4"]

    else:
        list_dir = [sys.argv[2]]
    for li in list_dir:
        thisChrom = li
        infile = thisChrom + ".fa"
        if thisChrom == 'chrY':
            infile = "chrY.psr.fa"
        chrom = []
        x = ""

        with open(infile, 'r') as INFILE:
            INFILE.readline()
            for x in INFILE:
                chrom.append(x.rstrip())

            x = "".join(chrom)

        chrom = None
        x = x.upper()

        for i in range(len(x) - read_length + 1):
            print(">" + li + "." + str(i))
            print(x[i:(i + read_length)])


if __name__ == "__main__":
    main()