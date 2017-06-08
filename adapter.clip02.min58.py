#!/usr/bin/env python

import sys
import time
import random
import math
import re


def main():

	adapter_seq = sys.argv[1]

	p = re.compile(adapter_seq)
	
	counter = 0
	prevRead = []
	
	for x in sys.stdin:
		aline = x.rstrip()
		
		if counter % 4 == 0:
			if len(prevRead) == 0:
				pass
			else:
				if readLength > 57:
					for pr in prevRead:
						print pr
				
			prevRead = []
			prevRead.append(aline)
		
		elif counter % 4 == 1:
			
			readLength = len(aline)
			m = p.search(aline)
			if m > 0:
				readLength = m.start()
			else:
				for i in range(4, 14):
					if aline[-i:] == adapter_seq[0:i]:
						readLength = readLength - i
			prevRead.append(aline[0:readLength])
			
		elif counter % 4 == 2:

			prevRead.append(aline)

		elif counter % 4 == 3:

			prevRead.append(aline[0:readLength])

		else:
			
			print "ERROR IN ADAPTER CLIP"
			
		counter += 1
			
	
if __name__ == "__main__":
	main()
