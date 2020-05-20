#!/usr/bin/env python3

import csv
import matplotlib.pyplot as plt
import numpy as np

from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform


if __name__ == "__main__":
	
	# with open("baumannii_dist.tsv", "r") as tsv_file:
	with open("sagalactiae_dist.tsv", "r") as tsv_file:

		read_tsv = csv.reader(tsv_file, delimiter="\t")

		dist_array = []
		for row in read_tsv:
			dist_array.append(row)
		
		tsv_file.close()

	dist_array = np.array(dist_array)
	label_array = dist_array[0,1:]

	dist_array = np.delete(dist_array, 0, axis=1)
	dist_array = np.delete(dist_array, 0, axis=0)

	dist_array = dist_array.astype(np.float)

	dists = squareform(dist_array)
	linkage_matrix = linkage(dists, "complete")
	dendrogram(linkage_matrix, labels=label_array, color_threshold=0.006)
	
	plt.title("Complete Link")
	plt.show()
