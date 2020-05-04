#!/usr/bin/env python

def process_file(file_name):
   
   lines=[]
   
   for line in open(file_name+".txt", "r"):
   
      if line.startswith("Gene ID") or line.strip() == '':
         continue
      lines.append(line.strip())
      
   return sorted(lines)

def main():
   sorted_lara=process_file("Laragenes")
   sorted_fc=process_file("FC_ensemblgenes")
   
   common_genes=list(set(sorted_lara).intersection(sorted_fc))
   for gene in common_genes:
      print(gene)

   common_gene_number=len(common_genes)
   print("common genes: " + str(common_gene_number))
   
   unique_gene_number_lara=len(sorted_lara)-common_gene_number
   print("unique genes in lara: "+ str(unique_gene_number_lara))
   
   unique_gene_number_fc=len(sorted_fc)-common_gene_number
   print("unique genes in FC_ensemblgenes: "+ str(unique_gene_number_fc))

if __name__  == "__main__":
   main()
