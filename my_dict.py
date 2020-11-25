#!/usr/bin/env python3

"""
Usage:
my_dict.py file <file_name>
my_dict.py all
my_dict.py add_to_all <file_name>

Options:
file <file_name>         Search for only the given text file given the following format. word: meaning \n word: meaning , etc.
all                      Search for the file added the words until now.
add_to_all <file_name>   Add text file to the whole set.
"""

from docopt import docopt
import os
import random


def get_dictionary():
	with open(file_name, "r") as file:

		for line in file.readlines():
			line = line.strip()
			word_meaning = line.split(":")
				
			if len(word_meaning) != 2:
				raise Exception("meaning is not given")
			else:
				dict_dict[word_meaning[0]] = word_meaning[1]

		file.close()

def check_all_file_existence():
	try:
		file = open(all_file_name)
	except FileNotFoundError:
		print("{} is not found. It is created.".format(all_file_name))
		with open(all_file_name, "w"):
			pass

def is_file_empty():
	return os.stat(file_name).st_size != 0

def copy_to_all_file():
	with open(file_name, "r") as in_file:
		check_all_file_existence()

		with open(all_file_name, "a") as out_file:
			for line in in_file.readlines():
				out_file.write(line)
			out_file.close()

		in_file.close()

def on_press(key):
	# print('{} is pressed'.format(key))
	if hasattr(key,'char'):
		if key.char == 'q':
			return False

def give_me_a_word():
	text = "Press any character to continue except q, q is for the quit:"

	while True:
		word, meaning = random.choice(list(dict_dict.items()))
		print("* What is the meaning of {}?\n".format(word))
		if input() == 'q':
			break
		print("{}\n".format(meaning))
		if input() == 'q':
			break
	

if __name__ == "__main__":

	args = docopt(__doc__)

	all_file_name = "all.txt"

	dict_dict = {}

	if args["file"]:
		file_name = args["<file_name>"]

	elif args["all"]:
		check_all_file_existence()

	elif args["add_to_all"]:
		file_name = args["<file_name>"]
		copy_to_all_file()

	if is_file_empty():
		get_dictionary()

	else:
		print(("{} is empty.").format(file_name))
		exit()

	if args["all"] or args["file"]:
		give_me_a_word()






