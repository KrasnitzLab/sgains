#!/usr/bin/env python

import sys
import time
from operator import itemgetter


def main():

	outfilename = "hg19.chrom.sizes.txt"

	chroms = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "M"]

	INFILE = 0
	OUTFILE = open(outfilename, "w")

	abspos = 0
	for x in chroms:
		infilename = "/filepath/chr" + x + ".fa"
		print infilename
		if INFILE:
			INFILE.close()
		INFILE = open(infilename, "r")
		INFILE.readline()
		chromlength = 0
		for y in INFILE:
			aline = y.rstrip()
			chromlength += len(aline)
		OUTFILE.write("chr" + x)
		OUTFILE.write("\t")
		OUTFILE.write(str(chromlength))
		OUTFILE.write("\t")
		OUTFILE.write(str(abspos))
		OUTFILE.write("\n")
		abspos += chromlength
		
	INFILE.close()
	OUTFILE.close()


if __name__ == "__main__":
	main()
