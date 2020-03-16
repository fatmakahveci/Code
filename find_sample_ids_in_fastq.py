#!/usr/bin/env python3.7.4

"""
Date: March 16, 2020
Author: Fatma Kahveci

Aim: To find the sample IDs using headers in FASTQ files

Usage:
	find_sample_ids_in_fastq.py --input=<fastq_file>

"""


## LIBRARIES ##

from docopt import docopt


## MAIN METHOD ##

if __name__ == "__main__":

	args = docopt(__doc__, version="0.6.2")
	
	sample_set = set()

	with open(args["--input"], 'r') as in_file:

		for line in in_file.readlines():

			if line.startswith('@') and 'length' in line:
				sample_set.add(line[1:].split('.')[0])
		
		in_file.close()

	print("Samples in FASTQ file = "+", ".join(sample_set))
