#!/usr/bin/env python3.7

from string import ascii_uppercase

if __name__ == "__main__":
	
	char_list = []

	for i in ascii_uppercase:
		char_list.append(i)
		open(i, 'a').close()

	with open("nucl_gb.accession2taxid", "r") as file:

		while True:
			
			line = file.readline()
			
			if not line:
				break
			
			if line[0] in char_list:
				temp_file = open(line[0],'a')
				temp_file.write(line)
				temp_file.close()

		file.close()
