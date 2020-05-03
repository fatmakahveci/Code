#!/usr/bin/env python3.7

from Bio import SeqIO
import sys

if __name__ == "__main__":

	input_file = str(sys.argv[1])
	word_search_in_headers = str(sys.argv[2])
	
	print("Usage: python search_in_fasta.py \"$input_file$\" \"$word$\" ")
	print("$"+word_search_in_headers + "$ is being searched...")
	
	out_file = open(input_file.split(".")[0]+".out."+input_file.split(".")[1], "w")

	for record in SeqIO.parse(input_file, "fasta"):
		id = str(record.id)
		seq = str(record.seq)
		if word_search_in_headers in record.description:
			out_file.write(">"+str(record.description)+"\n"+str(record.seq)+"\n")

	out_file.close()
