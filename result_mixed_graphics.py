#!/usr/bin/env python3.8.2


import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


if __name__ == "__main__":

	(f, (ax1,ax2)) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [2.5, 1]})

	## FIRST PLOT ##

	category_names = ['RA(%70) - First sample', 'RA(%20) - Second sample', 'RA(%10) - Third sample']

	results = {
		'BWA': [66.52, 20.10, 13.38],
		'SBG': [66.80, 20.02, 13.18],
		'VG': [66.48, 19.98, 13.54]
	}

	actual_results = {
		'BWA': [70.0, 90.0, 100.0],
		'SBG': [70.0, 90.0, 100.0],
		'VG': [70.0, 90.0, 100.0]
	}

	labels = list(results.keys())
	actual_results_list = np.array(list(actual_results.values()))
	data = np.array(list(results.values()))
	data_cum = data.cumsum(axis=1)
	category_colors = plt.get_cmap('RdYlBu')(np.linspace(0.15, 0.85, data.shape[1]))

	ax1.invert_yaxis()
	ax1.xaxis.set_visible(False)
	ax1.set_xlim(0, np.sum(data, axis=1).max())

	for i, (colname, color) in enumerate(zip(category_names, category_colors)):
		widths = data[:, i]
		actual_value = actual_results_list[:,i]
		starts = data_cum[:, i] - widths
		ax1.barh(labels, widths, left=starts, height=0.5, label=colname, color=color)

		xcenters = starts + widths / 2

		r, g, b, _ = color
		# text_color = 'white' if r * g * b < 0.5 else 'darkgrey'
		text_color = 'black'
		for y, (x, c) in enumerate(zip(xcenters, widths)):
			ax1.text(x, y, str(float(c)), ha='center', va='center', color=text_color)

	# ax1.legend(ncol=len(category_names), bbox_to_anchor=(0, 1), loc='lower left', fontsize='small')

	ax1.axvline(70, color="black", ymin=0.05, ymax=0.225, linestyle='dotted')
	ax1.axvline(70, color="black", ymin=0.41, ymax=0.59, linestyle='dotted')
	ax1.axvline(70, color="black", ymin=0.78, ymax=0.96, linestyle='dotted')

	ax1.axvline(90, color="black", ymin=0.05, ymax=0.225, linestyle='dotted')
	ax1.axvline(90, color="black", ymin=0.41, ymax=0.59, linestyle='dotted')
	ax1.axvline(90, color="black", ymin=0.78, ymax=0.96, linestyle='dotted')

	left, width = .25, .5
	bottom, height = .25, .5
	right = left + width
	top = bottom + height

	ax1.text(0.35*(left+right), 0.72*(bottom+top), 'RA:70%',
		horizontalalignment='center',
		verticalalignment='center',
		fontsize=8, color='black',
		transform=ax1.transAxes)

	ax1.text(0.80*(left+right), 0.72*(bottom+top), 'RA:20%',
		horizontalalignment='center',
		verticalalignment='center',
		fontsize=8, color='black',
		transform=ax1.transAxes)

	ax1.text(0.95*(left+right), 0.72*(bottom+top), 'RA:10%',
		horizontalalignment='center',
		verticalalignment='center',
		fontsize=8, color='black',
		transform=ax1.transAxes)

	ax1.text(0.35*(left+right), 0.36*(bottom+top), 'RA:70%',
		horizontalalignment='center',
		verticalalignment='center',
		fontsize=8, color='black',
		transform=ax1.transAxes)

	ax1.text(0.80*(left+right), 0.36*(bottom+top), 'RA:20%',
		horizontalalignment='center',
		verticalalignment='center',
		fontsize=8, color='black',
		transform=ax1.transAxes)

	ax1.text(0.95*(left+right), 0.36*(bottom+top), 'RA:10%',
		horizontalalignment='center',
		verticalalignment='center',
		fontsize=8, color='black',
		transform=ax1.transAxes)

	ax1.text(0.35*(left+right), 0.02*(bottom+top), 'RA:70%',
		horizontalalignment='center',
		verticalalignment='center',
		fontsize=8, color='black',
		transform=ax1.transAxes)

	ax1.text(0.80*(left+right), 0.02*(bottom+top), 'RA:20%',
		horizontalalignment='center',
		verticalalignment='center',
		fontsize=8, color='black',
		transform=ax1.transAxes)

	ax1.text(0.95*(left+right), 0.02*(bottom+top), 'RA:10%',
		horizontalalignment='center',
		verticalalignment='center',
		fontsize=8, color='black',
		transform=ax1.transAxes)

	ax1.grid(which='major', linestyle='-', linewidth='0.05', alpha=0.3, color='black')
	ax1.grid(which='minor', linestyle='--', linewidth='0.001', alpha=0.1, color='gray')
	ax1.grid(True)

	## SECOND PLOT ##
	n_clusters = 5
	# sample = np.arange(0, 6)
	
	# bwa_sa_before_est = np.array([40, 66, 57, 31, 32])
	# bwa_sa_after_est = np.array([94, 127, 121, 83, 80])
	
	# sbg_sa_before_est = np.array([26, 32, 32, 21, 31])
	# sbg_sa_after_est = np.array([216, 284, 346, 384, 333])

	# vg_sa_before_est = np.array([227, 239, 234, 220, 221])
	# vg_sa_after_est = np.array([280, 300, 297, 271, 266])

	index = np.arange(1, 2)
	
	bar_width = 0.2

	ax2.set(xlim=(0.5, 1.5), ylim=(0, 700))

	ax2.bar(index - bar_width, 45.2, bar_width, color='navy', edgecolor='black', label='BWA')
	ax2.bar(index - bar_width, 101.0, bar_width, color='dodgerblue', bottom = 45.2, edgecolor='black', label='BWA + Estimation')
	
	ax2.bar(index, 28.4, bar_width, color='darkgreen', edgecolor='black', label='SBG')
	ax2.bar(index, 312.6, bar_width, color='seagreen', bottom = 28.4, edgecolor='black', label='SBG + Estimation')
	
	ax2.bar(index + bar_width, 228.2, bar_width, color='darkorange', edgecolor='black', label='VG')
	ax2.bar(index + bar_width, 282.8, bar_width, color='orange', bottom = 228.2, edgecolor='black', label='VG + Estimation')
	
	# ax2.title.set_text("Run time")
	ax2.set_ylabel('Run time (sec)')
	ax2.set_xticks(np.arange(0.8, 1.4, 0.2))
	ax2.set_yticks(np.arange(0, 700, 50))
	ax2.set_xticklabels(['BWA', 'SBG', 'VG'], rotation=90)
	ax2.minorticks_on()
	ax2.grid(which='major', linestyle='-', linewidth='0.05', alpha=0.3, color='black')
	ax2.grid(which='minor', linestyle='--', linewidth='0.001', alpha=0.1, color='gray')
	ax2.grid(True)

	ax2.legend(loc='upper right', fontsize='x-small')
	
	plt.tight_layout()
	# plt.savefig('mixed_10_3_2.eps', format='eps')
	plt.show()
