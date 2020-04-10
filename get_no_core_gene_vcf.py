#!/usr/bin/env python3.7


if __name__ == "__main__":
	
	n = 10

	in_vcf_file = open("strain_graph.vcf", "r")
	out_vcf_file = open("strain_graph"+str(n)+".vcf", "w")

	in_fasta_file = open("baumannii.fasta", "r")
	out_fasta_file = open("baumannii"+str(n)+".fasta", "w")

	contig_idx = 0

	added_gene_list = list()

	for line in in_vcf_file.readlines():

		if line.startswith("##contig"):
	
			if contig_idx < n:
	
				out_vcf_file.write(line)
				contig_idx += 1
	
		elif line.startswith("#"):
			out_vcf_file.write(line)

		else:
			
			gene_id = line.strip("\n").split("\t")[0]

			if gene_id not in added_gene_list:
				added_gene_list.append(gene_id)

			if len(added_gene_list) <= n:
				out_vcf_file.write(line)

			elif len(added_gene_list) == n+1:
				break

	vcf_gene_no = int(gene_id.split("ACIBA")[1])

	for line in in_fasta_file.readlines():

		if line.startswith(">"):

			fasta_gene_no = int(line.split("\t")[0].split("ACIBA")[1])

		if fasta_gene_no < vcf_gene_no:
			out_fasta_file.write(line)

	in_vcf_file.close()
	out_vcf_file.close()

	in_fasta_file.close()
	out_fasta_file.close()
