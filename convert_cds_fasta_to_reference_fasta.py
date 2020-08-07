#!/usr/bin/env python3.7

"""
Aim: It converts coding sequences' fasta file into reference fasta file.
Requirements: pip3 install Bio, docopt

Usage:
    convert_cds_into_ref_fasta.py --cds_fasta <cds_fasta_file_name> --output_fasta <output_fasta_file_name>
    convert_cds_into_ref_fasta.py --help

Options:
    -h, --help                                  Show this help message and exit.
    --cds_fasta <cds_fasta_file_name>           VCF file containing snvs and strain existence info [default: ].
	--output_fasta <output_fasta_file_name>		Output file name [default: ].
"""

from Bio import SeqIO
from docopt import docopt


if __name__ == "__main__":

	args = docopt(__doc__, version="0.6.2")
	
	in_file = args["--cds_fasta"]
	out_file = args["--output_fasta"]

	fasta_list = []

	with open(in_file, 'r') as in_file:

		for record in SeqIO.parse(in_file, "fasta"):
			fasta_list.append(str(record.seq))

	print('N'.join(fasta_list))
