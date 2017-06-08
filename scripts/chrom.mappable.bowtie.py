#!/usr/bin/env python

import re
import sys
import random
import string

def main():

	INFILE = open("mappable.regions.sorted.txt", "r")
	OUTFILE = open("chrom.mappable.bowtie.txt", "w")

	prevChrom = "chr1"
	prevChromMappable = 0

	OUTFILE.write("chrom")
	OUTFILE.write("\t")
	OUTFILE.write("num_mappable_pos")
	OUTFILE.write("\n")

	for x in INFILE:
		arow = x.rstrip().split("\t")
		thisChrom = arow[0]
		thisLength = int(arow[2]) - int(arow[1])

		if thisChrom == prevChrom:
			prevChromMappable += thisLength
		else:
			OUTFILE.write(prevChrom)
			OUTFILE.write("\t")
			OUTFILE.write(str(prevChromMappable))
			OUTFILE.write("\n")
			prevChromMappable = thisLength
		prevChrom = thisChrom

	OUTFILE.write(prevChrom)
	OUTFILE.write("\t")
	OUTFILE.write(str(prevChromMappable))
	OUTFILE.write("\n")
	
	INFILE.close()	
	OUTFILE.close()


if __name__ == "__main__":
	main()
