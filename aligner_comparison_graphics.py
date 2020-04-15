#!/usr/bin/env python3.8.2


import matplotlib.pyplot as plt
import numpy as np
import math


if __name__ == "__main__":

	n_clusters = 3

	f, ((ax1, ax2)) = plt.subplots(1, 2)

	## FIRST PLOT ##

	wg_throughput = [19549, 40796, 5605]
	cg_throughput = [6173, 16035, 492]

	index = np.arange(1,n_clusters+1)
	bar_width = 0.1
	opacity = 1.0

	ax1.bar(index - bar_width, wg_throughput, 2 * bar_width, alpha=opacity, label='wg')
	ax1.bar(index + bar_width, cg_throughput, 2 * bar_width, alpha=opacity, label='cg')

	ax1.title.set_text("Throughput")
	ax1.set_ylabel('Number of mapped reads per second')
	ax1.set_xticks(np.arange(1, 4))
	ax1.set_xticklabels(['BWA', 'SBG', 'VG'])
	ax1.minorticks_on()
	ax1.grid(which='major', linestyle='-', linewidth='0.05', alpha=0.3, color='black')
	ax1.grid(which='minor', linestyle='--', linewidth='0.001', alpha=0.1, color='gray')
	ax1.grid(True)
	ax1.legend(loc='upper left')

	## SECOND PLOT ##

	wg_coverage = [503, 503, 504]
	cg_coverage = [538, 498, 682]

	ax2.bar(index - bar_width, wg_coverage, 2 * bar_width, alpha=opacity, label='wg')
	ax2.bar(index + bar_width, cg_coverage, 2 * bar_width, alpha=opacity, label='cg')

	ax2.title.set_text("Mappability")
	ax2.set_ylabel('Coverage')
	ax2.set_xticks(np.arange(1, 4))
	ax2.set_xticklabels(['BWA', 'SBG', 'VG'])
	ax2.minorticks_on()
	ax2.grid(which='major', linestyle='-', linewidth='0.05', alpha=0.3, color='black')
	ax2.grid(which='minor', linestyle='--', linewidth='0.001', alpha=0.1, color='gray')
	ax2.grid(True)
	ax2.legend(loc='upper left')

	plt.tight_layout()
	plt.savefig('aligner_comparison.pdf', format='pdf')
	
