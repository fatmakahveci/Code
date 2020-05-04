#!/bin/env/python

def main():
	for sample in ['ERR249240','ERR249291','SRR1200553','SRR1314629','SRR3131131','SRR3146384']:
		file_content=open('/Users/fatma/Desktop/bacgentrack/data/wgs/'+sample+'.sample_specific.bed','r').read().split('\n')
		if file_content != '':
			second_col=[]
			third_col=[]

			for line in file_content:
				fields=line.split('\t')
				second_col.append(int(fields[1]))
				third_col.append(int(fields[2]))
			print str(sum(third_col)-sum(second_col))

if __name__=="__main__":
	main()
