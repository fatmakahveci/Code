#!/bin/python

import sys

directory="/home/fatma.balci/projects/gibbon/delly/output/del/"

tp_file_name="del_gene_t_p.bed"
nt_file_name="del_gene_n_t.bed"

output_file_name="delly_del_n_p_t_result.bed"

try:
    nt_file=open(directory+nt_file_name,"r")
    tp_file=open(directory+tp_file_name,"r")

    output_file = open(directory + output_file_name, "w")
except(OSError, IOError) as e:
    sys.exit(e)

tp_line_number=0;

tp_file_content=tp_file.readlines(); nt_file_content=nt_file.readlines()

tp_total_line_number=len(tp_file_content); nt_total_line_number=len(nt_file_content)
output_file.write("chr \t start \t end \t GT(n) \t GT(t) \t GT(p) \t geneChr \t geneStart \t geneEnd \t geneName \t %n \t %t \t %p \n")

while tp_line_number < tp_total_line_number:

    tpfields = tp_file_content[tp_line_number].strip().split('\t')

    nt_line_number=0

    while nt_line_number < nt_total_line_number:
        ntfields = nt_file_content[nt_line_number].strip().split('\t')
        if tpfields[3] == ntfields[3] and tpfields[7] == ntfields[7] and tpfields[8] == ntfields[8] and \
                        tpfields[9] == ntfields[9]:
            output_file.write(ntfields[6] + '\t' + ntfields[7] + '\t' + ntfields[8] + '\t' + ntfields[9] + '\t' + \
                              ntfields[10] + '\t' + tpfields[10] + '\t' + ntfields[0] + '\t' + ntfields[1] + '\t' + \
                              ntfields[2] + '\t' + ntfields[3] + '\t' + ntfields[4] + '\t' + ntfields[5] + '\t' + tpfields[5] + '\n')
        nt_line_number+=1
    tp_line_number+=1


nt_file.close()
tp_file.close()

output_file.close()
