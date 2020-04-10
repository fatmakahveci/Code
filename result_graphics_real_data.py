#!/usr/bin/env python3.7

import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
	n_clusters = 5
		
	bwa_aligner = []
	bwa_all = []
	for values in {"c7":[1, 883, 230, 574, 2291], "c8":[1, 2546, 87, 141, 652], "c9":[1, 472, 132, 143, 255], "c10":[1, 443, 82, 88, 808], "c12":[1, 945, 252, 370, 1049]}.values():
		bwa_aligner.append(round(sum(values[:-1]) / 60.0, 2))
		bwa_all.append(round(sum(values) / 60.0, 2))

	sbg_aligner = []
	sbg_all = []
	for values in {"c7":[351, 916, 2327], "c8":[468, 289, 616], "c9":[280, 145, 249], "c10":[152, 115, 794], "c12":[110, 193, 1024]}.values():
		sbg_aligner.append(round(sum(values[:-1]) / 60.0, 2))
		sbg_all.append(round(sum(values) / 60.0, 2))

	vg_aligner = []
	vg_all = []
	for values in {"c7":[1, 43, 19035, 292, 2406], "c8":[1, 43, 6801, 117, 695], "c9":[1, 43, 2800, 72, 276], "c10":[1, 43, 5339, 85, 916], "c12":[1, 43, 10548, 184, 1056]}.values():
		vg_aligner.append(round(sum(values[:-1]) / 60.0, 2))
		vg_all.append(round(sum(values) / 60.0, 2))

	index = np.arange(n_clusters)
	bar_width = 0.3
	opacity = 0.8

	f, ((ax1, ax2)) = plt.subplots(1, 2, sharex='col', sharey='row')

	ax1.bar(index, sbg_aligner, bar_width, alpha=opacity, color='b', label='bwa')
	ax1.bar(index + bar_width, bwa_aligner, bar_width, alpha=opacity, color='g', label='SBG')
	ax1.bar(index + 2 * bar_width, vg_aligner, bar_width, alpha=opacity, color='y', label='VG')
	ax1.title.set_text("Aligner")
	ax1.set_xlabel('cluster ID')
	ax1.set_ylabel('run time in minutes')
	ax1.grid(which='major', linestyle='-', linewidth='0.05', color = "black")
	ax1.grid(which='minor', linestyle='--', linewidth='0.01', color = "black")
	ax1.grid(True)
	ax1.legend()

	ax2.bar(index, sbg_all, bar_width, alpha=opacity, color='b', label='bwa')
	ax2.bar(index + bar_width, bwa_all, bar_width, alpha=opacity, color='g', label='SBG')
	ax2.bar(index + 2 * bar_width, vg_all, bar_width, alpha=opacity, color='y', label='VG')
	ax2.title.set_text("Aligner and Estimation")
	ax2.set_xlabel('cluster ID')
	ax2.set_ylabel('run time in minutes')
	ax2.grid(which='major', linestyle='-', linewidth='0.05', color = "black")
	ax2.grid(which='minor', linestyle='--', linewidth='0.01', color = "black")
	ax2.grid(True)
	ax2.legend()

	plt.tight_layout()
	plt.show()
	# plt.savefig('common_labels.png', dpi=300)
