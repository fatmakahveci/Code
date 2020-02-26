#!/usr/bin/env python3.7

# usage: python3.7 add_snps_per_given_number_of_base.py --fasta ${file_name}.fasta --interval ${per_given_number_of_bases}

import argparse
from Bio import SeqIO
from random import seed, randint


def change_base(base):

	base_dict = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
	
	return str(base_dict[base])


if __name__ == "__main__":

	parser = argparse.ArgumentParser()

	parser.add_argument("--fasta", help = "Input Fasta File")
	parser.add_argument("--interval", help = "SNP per given interval")

	args = parser.parse_args()

	interval = int(args.interval)

	out_file = open(args.fasta.strip(".fasta")+"_"+str(interval)+".fasta","w")
	
	for fasta_seq in SeqIO.parse(args.fasta, "fasta"):

		seq = str(fasta_seq.seq)
		seq_id = str(fasta_seq.id)

		seq_len = len(seq)
		
		last_idx = seq_len

		while last_idx >= 0:

			random_num = randint(0, interval)

			if seq_len > interval:
				for i in range(50):
					random_num = randint(0, interval)
					seq = seq[:random_num] + change_base(seq[random_num]) + seq[random_num+1:]

			else:
				for i in range(seq_len % 20):
					random_num = randint(0, interval)
				random_num = randint(0, seq_len)

			seq = seq[:random_num] + change_base(seq[random_num]) + seq[random_num+1:]

			last_idx -= interval

		out_file.write(">"+seq_id+"\n"+seq+"\n")

	out_file.close()
