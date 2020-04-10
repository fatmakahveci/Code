#!/usr/bin/env python3.7


if __name__ == "__main__":
	
	n = 10

	in_vcf_file = open("strain_graph.vcf", "r")
	out_vcf_file = open("strain_graph"+str(n)+".vcf", "w")

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
			
			out_vcf_file.write(line)

			if gene_id not in added_gene_list:
				added_gene_list.append(gene_id)

			if len(added_gene_list) == n:
				break
	
	in_vcf_file.close()
	out_vcf_file.close()
