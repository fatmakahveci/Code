#!/usr/bin/env python3.7

import matplotlib.pyplot as plt
import numpy as np
import math

if __name__ == "__main__":
	n_clusters = 5
	
	clusters_all = [x/1000000 for x in map(float, [17876054, 2966490, 921602, 4665824, 9269392])]

	bwa_mapped = [x/1000000 for x in map(float, [5792218, 1015127, 345921, 1675525, 2400259])]
	sbg_mapped = [x/1000000 for x in map(float, [5556580, 884753, 311250, 1566679, 2267545])]
	vg_mapped = [x/1000000 for x in map(float, [7380929, 1664377, 546170, 2285533, 3306067])]
	
	f, ((ax1, ax2)) = plt.subplots(1, 2)

	index = np.arange(1,n_clusters+1)
	bar_width = 0.1
	opacity = 0.3
	ax1.bar(index - bar_width, clusters_all, bar_width, alpha=opacity, color='b', label='all')
	opacity = 0.8
	ax1.bar(index - bar_width, bwa_mapped, bar_width, alpha=opacity, color='b', label='SBG')
	opacity = 0.3
	ax1.bar(index, clusters_all, bar_width, alpha=opacity, color='g')
	opacity = 0.8
	ax1.bar(index, sbg_mapped, bar_width, alpha=opacity, color='g', label='bwa')
	opacity = 0.3
	ax1.bar(index + bar_width, clusters_all, bar_width, alpha=opacity, color='y')
	opacity = 0.8
	ax1.bar(index + bar_width, vg_mapped, bar_width, alpha=opacity, color='y', label='VG')

	ax1.title.set_text("Read Mapping")
	ax1.set_xlabel('cluster ID')
	ax1.set_ylabel('number of mapped reads (x1,000,000)')
	ax1.set_xticks(np.arange(1,6,1))
	ax1.grid(which='major', linestyle='-', linewidth='0.05', color = "black")
	ax1.grid(which='minor', linestyle='--', linewidth='0.01', color = "black")
	ax1.grid(True)
	ax1.legend(loc='upper right')

	genome_size = [1800, 889.8, 259.8, 699.4, 936.2]
	index = np.arange(1,n_clusters+1)
	bar_width = 0.1
	opacity = 0.8
	ax2.bar(index - bar_width, genome_size, bar_width, alpha=opacity, color='b')

	ax2.title.set_text("Genome size")
	ax2.set_xlabel('cluster ID')
	ax2.set_ylabel('genome size (Mbp)')
	ax2.set_xticks(np.arange(1,6,1))
	ax2.grid(which='major', linestyle='-', linewidth='0.05', color = "black")
	ax2.grid(which='minor', linestyle='--', linewidth='0.01', color = "black")
	ax2.grid(True)
	ax2.legend(loc='upper right')

	plt.tight_layout()
	plt.show()
