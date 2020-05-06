#!/usr/bin/env python3

import os, string


if __name__=="__main__":

	upper_list = list(string.ascii_uppercase)

	for letter in upper_list:
		cmd = "touch "+letter
		os.system(cmd)

	for file in ["acc_tax_gb.txt", "acc_tax.txt"]:
	
		with open(file, 'r') as in_file:

			while True:
			
				line = in_file.readline().strip()
				cmd = 'echo \"'+line+'\" >> '+ line[0]
				os.system(cmd)

				if not line:
					break

		in_file.close()
