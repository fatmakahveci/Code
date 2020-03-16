#!/usr/bin/env python3.7

"""
Date: March 16, 2020
Author: Fatma Kahveci

Aim: This code aims to add info and extract the columns belonging to the given samples from species' reference vcf file.

Usage:
	correct_vcf.py --vcf=<input_file_name> --out=<output_file_name> [--correct_info] [--remove_same_gt_line] (--sample=<sample> | --samples_in_fastq=<fastq_file_name>)
	correct_vcf.py --help

Options:
	-h --help					Show this help message and exit.
	--vcf=input_file_name		Input VCF file name
	--out=output_file_name		Output VCF file name
	--sample=sample_list		List of samples to be extracted. [default: None]
	--correct_info				Is info field will be corrected? [default: False].
	--remove_same_gt_line		Is the line with same gt for each sample removed? [default: False].
	--samples_in_fastq			Take sample names from FASTQ file. [default: None]
"""


## LIBRARIES ##

from docopt import docopt
import sys

## METHODS ##

def correct_line(line): # add required information instead of '.' and merge the fields '\t'

	fields = get_fields(line)

	if args["--correct_info"]:
		fields[7] = "AC=10;AN=60;"
		fields[8] = "GT"

	return '\t'.join(fields)

def get_fields(line):
	return line.strip('\n').split('\t')

def get_gt_fields(line):
	return get_fields(line)[no_info_field:]

def get_info_fields(line):
	return get_fields(line)[:no_info_field]

def get_sample_col_idx(line, sample_list): # index of given sample in #CHROM header
	
	idx_list = list()

	is_item_in_list = lambda item, item_list: item_list.index(item) if item in item_list else None
		
	for sample in sample_list:

		idx = is_item_in_list(sample, get_fields(line))
			
		if not idx is None:
			idx_list.append(idx)

	if len(idx_list) == 0:

		print(",".join(sample_list) + " column do not exist in vcf file.")
		
		sys.exit(0)

	return idx_list

def is_header(line): # check whether if the line is header

	if line.startswith('##'):
		return True
	
	return False

def is_info_header(line): # Is the line starts with #CHROM

	if line.startswith("#CHROM"):
		return True

	return False

def is_sample_genome_line(line): # check whether if ##SAMPLE line among header lines
	if line.startswith('##SAMPLE'):
		return True
	
	return False

def extract_sample_from_sample_line(line): # ##SAMPLE=<ID=Strain_Cluster,GENOMES=SampleID1,SampleID2> -> ##SAMPLE=<ID=Strain_Cluster,GENOMES=SampleID1>
	
	for sample in sample_list:

		line = line.replace(sample,"").replace(";;",";").replace(";>",">")

	return line

def extract_sample_column(line): # #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SampleID1 SampleID2 -> #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SampleID2 
	
	fields = get_fields(line)
	
	for idx in sorted(sample_idx_list, reverse=True):
		try:
			del fields[idx]
		except:
			return ""

	return "\t".join(fields)

def get_sample_names_from_fastq():

	sample_set = set()

	with open(args["--samples_in_fastq"], 'r') as in_file:

		for line in in_file.readlines():

			if line.startswith('@') and 'length' in line:
				sample_set.add(line[1:].split('.')[0])
		
		in_file.close()

	return list(sample_set)

def extract_line_with_same_gt(line): # extract the lines that are incapable of representation

	gts = get_gt_fields(line)
	gt_to_compare = gts[0]

	for gt in gts:
		
		if gt != gt_to_compare:
			return line+'\n'
	
	return ""

def get_sample_list():

	if not args["--samples_in_fastq"] is None:
		return get_sample_names_from_fastq()

	else:
		return args["--sample"].split(",")



### MAIN METHOD ###

if __name__ == "__main__":
	
	args = docopt(__doc__, version="0.6.2")

	sample_list = get_sample_list()

	out_file = open(args["--out"], 'w')

	no_info_field = 9

	with open(args["--vcf"], 'r') as in_file:

		for line in in_file.readlines():

			if line == "":
				continue
			
			elif is_sample_genome_line(line): # ##SAMPLE
				out_file.write(extract_sample_from_sample_line(line))

			elif is_header(line): # All header files except ##SAMPLE and #CHROM
				out_file.write(line)

			elif is_info_header(line): # #CHROM HEADER
				
				sample_idx_list = get_sample_col_idx(line, sample_list)
				out_file.write(extract_sample_column(line)+'\n')

			else: # remove sample column

				if (args["--remove_same_gt_line"] and extract_line_with_same_gt(correct_line(line)) != "") or not args["--remove_same_gt_line"]:
					
					new_line = extract_sample_column(line)
					
					if new_line != "":
						out_file.write(new_line+'\n')


		in_file.close()

	out_file.close()

	
