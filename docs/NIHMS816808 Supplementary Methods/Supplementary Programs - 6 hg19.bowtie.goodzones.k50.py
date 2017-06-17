#!/usr/bin/env python

import re
import sys
import random
import string

def main():

	OUTFILE = open("/filepath/hg19.goodzones.bowtie.k50.bed", "w")
	
	inGood = False
	prevChrom = ""
	thisStart = -1
	thisEnd = -1
	for i in range(21):
		print i
		sys.stdout.flush()
		INFILE = open("/filepath/sequence.part." + str(i) + ".k50.sam", "r")
		for x in INFILE:
			if x[0] == "@":
				continue
			arow = x.rstrip().split("\t")
			thisId = arow[0]
			thisIdChrom = thisId.split(".")[0]
			thisIdChrompos = int(thisId.split(".")[1])
			thisCode = arow[1]
			thisChrom = arow[2]
			thisChrompos = int(arow[3])
			if thisChrom == prevChrom:
				pass
			else:
				if inGood:
					OUTFILE.write(prevChrom)
					OUTFILE.write("\t")
					OUTFILE.write(str(thisStart))
					OUTFILE.write("\t")
					OUTFILE.write(str(thisEnd))
					OUTFILE.write("\n")
					inGood = False

			if thisCode == "0" and (thisChrompos - 1) == thisIdChrompos:
				if inGood:
					thisEnd = thisChrompos
				else:
					inGood = True
					thisStart = thisIdChrompos
					thisEnd = thisChrompos
			else:
				if inGood:
					OUTFILE.write(thisChrom)
					OUTFILE.write("\t")
					OUTFILE.write(str(thisStart))
					OUTFILE.write("\t")
					OUTFILE.write(str(thisEnd))
					OUTFILE.write("\n")
					inGood = False
			prevChrom = thisChrom

		INFILE.close()
	
	OUTFILE.close()


if __name__ == "__main__":
	main()
