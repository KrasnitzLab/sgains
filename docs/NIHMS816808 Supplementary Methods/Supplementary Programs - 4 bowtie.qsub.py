#!/usr/bin/env python

import sys
import time
from operator import itemgetter
import os

def main():

	bowtie_dir = "/filepath/bowtie-0.12.7"
	samtools_dir = "/filepath/samtools-0.1.16"
	project_dir = "/filepath/simulation.hg19.01"

	for i in range(21):
		thisSeqfile = project_dir + "/sequence.part." + str(i) + ".k50.txt"
		thisSamfile = project_dir + "/sequence.part." + str(i) + ".k50.sam"
		
                qsubFile = project_dir + "/bowtie.k50." + str(i) + ".qsub"
                QSUB = open(qsubFile, "w")
                outline = '#$ -S /bin/bash\n'
                QSUB.write(outline)
                outline = '#$ -l virtual_free=3.8G\n'
                QSUB.write(outline)
                outline = bowtie_dir + "/bowtie -S -t -n 2 -e 70 -m 1 --best --strata --solexa1.3-quals hg19 " + thisSeqfile + " " + thisSamfile + "\n"
                QSUB.write(outline)
                QSUB.close()

                time.sleep(1)
                thisCommand = "chmod 755 " + qsubFile
                os.system(thisCommand)
                thisCommand = "qsub " + qsubFile
                os.system(thisCommand)


if __name__ == "__main__":
	main()
