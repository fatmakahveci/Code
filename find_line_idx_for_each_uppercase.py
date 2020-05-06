#!/usr/bin/env python3

import string


if __name__=="__main__":

	upper_list = list(string.ascii_uppercase)

	upper_idx_dict = {}

	with open("acc_tax_gb.txt", "r") as f:

		for idx, line in enumerate(f, 1):
			
			first_char = line[0]

			if (first_char in upper_list) and not (first_char in upper_idx_dict):

				upper_idx_dict[first_char] = idx

		f.close()

	print(upper_idx_dict)
