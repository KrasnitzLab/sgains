#!/usr/bin/env python

import sys
import time
import string
import os

def main():

	infilename = sys.argv[1]
	outdir = sys.argv[2]
	barcodefilename = sys.argv[3]
	lane = sys.argv[4]

	bcLength = 7	
	
	BARCODE = open(barcodefilename, "r")
	INFILE = open(infilename, "r")
	
	barcode = dict()
	for x in BARCODE:
		arow = x.rstrip().split()
		thisBarcodeId = arow[0]
		thisBarcode = arow[1]
		if barcode.has_key(thisBarcode):
			print "ERROR: Duplicate barcode."
			return -1
		else:
			barcode[thisBarcode] = thisBarcodeId
	
	print barcode		

	OUTFILES = dict()
	for k, v in barcode.items():
		thisOutfilename = outdir + "/s_" + lane + "_sequence." + v + ".txt"
		OUTFILES[v] = open(thisOutfilename, "w")
		
	thisRead = None
	thisRead = []
	
	thisReadId = ""
	thisBcId = ""
	counter = 0
	
	for aline in INFILE:
		aline = aline.rstrip()

		if counter % 4 == 0:

			##  WRITE PREVIOUS READ  ##
			if thisBcId != "":
				for j in range(len(thisRead)):
					OUTFILES[thisBcId].write(thisRead[j])
					OUTFILES[thisBcId].write("\n")
					OUTFILES[thisBcId].flush()


			##  START NEW READ  ##
			thisBcId = ""
			thisRead = None
			thisRead = []
			thisRead.append(aline)
			

		elif counter % 4 == 1:

			thisBc = aline[0:bcLength]
			if barcode.has_key(thisBc):
				thisBcId = barcode[thisBc]
				thisRead.append(aline)


		elif counter % 4 == 2:
			thisRead.append(aline)

		else:
			thisRead.append(aline)

		counter += 1
		
	INFILE.close()
	for i in range(len(OUTFILES)):
		if OUTFILES[i]:
			OUTFILES[i].close()
	

if __name__=="__main__":
	main()
