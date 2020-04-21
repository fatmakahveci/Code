#!/usr/bin/env python3.8.2


import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


if __name__ == "__main__":

	(f, (ax1,ax2)) = plt.subplots(1, 2)

	## FIRST PLOT ##

	n_clusters = 5
	
	bwa_content1 = np.array([100, 100, 100, 100, 100])
	sbg_content1 = np.array([100, 100, 100, 100, 100])
	vg_content1 = np.array([100, 100, 100, 100, 100])

	bwa_content2 = np.array([100, 100, 100, 100, 100])
	sbg_content2 = np.array([100, 100, 100, 100, 100])
	vg_content2 = np.array([22.7, 24.6, 18.8, 9.5, 24.3])

	index = np.arange(1, n_clusters + 1)
	
	bar_width = 0.2

	ax1.set(xlim=(0.5, 5.5), ylim=(0, 100))

	ax1.bar(index - bar_width, bwa_content1, bar_width, color='navy', edgecolor='black')
	ax1.bar(index - bar_width, bwa_content2, bar_width, color='dodgerblue', bottom = bwa_content1, edgecolor='black')
	
	ax1.bar(index, sbg_content1, bar_width, color='darkgreen', edgecolor='black')
	ax1.bar(index, sbg_content2, bar_width, color='seagreen', bottom = sbg_content1, edgecolor='black')
	
	ax1.bar(index + bar_width, vg_content1, bar_width, color='darkorange', edgecolor='black')
	ax1.bar(index + bar_width, vg_content2, bar_width, color='orange', bottom = vg_content1, edgecolor='black')
	
	for i in range(1,6):
		ax1.text(i-0.2, 1.0, 'BWA', color='white', fontweight='bold',ha='center', va='bottom', fontsize='x-small', rotation=90)
		ax1.text(i, 9.0, 'SBG', color='white',fontweight='bold',ha='center', va='bottom', fontsize='x-small', rotation=90)
		ax1.text(i+0.2, 16.0, 'VG', color='white',fontweight='bold',ha='center', va='bottom', fontsize='x-small', rotation=90)

	# ax1.legend(bbox_to_anchor=(1.05, 1.0), fontsize='small', loc='upper left')
	# ax1.title.set_text("Strain Content\n(Single Sample)")
	# ax1.set_xlabel('Experiments')
	ax1.set_ylabel('Content (%)')
	ax1.set_xticks(np.arange(1, 6))
	ax1.set_yticks(np.arange(0, 110, 10))
	ax1.set_xticklabels(['S1', 'S2', 'S3', 'S4', 'S5'])
	ax1.minorticks_on()
	ax1.grid(which='major', linestyle='-', linewidth='0.05', alpha=0.3, color='black')
	ax1.grid(which='minor', linestyle='--', linewidth='0.001', alpha=0.1, color='gray')
	ax1.grid(True)

	## SECOND PLOT ##
	
	# sample = np.arange(0, 6)
	
	bwa_sa_before_est = np.array([1688, 2775, 748, 614, 1568])
	bwa_sa_after_est = np.array([3979, 3427, 1003, 1422, 2617])
	
	sbg_sa_before_est = np.array([1267, 757, 425, 267, 303])
	sbg_sa_after_est = np.array([3594, 1373, 674, 1061, 1327])
		
	vg_sa_before_est = np.array([19371, 6962, 2916, 5468, 10776])
	vg_sa_after_est = np.array([21777, 7657, 3192, 6384, 11832])

	index = np.arange(1, n_clusters + 1)
	
	bar_width = 0.2

	ax2.set(xlim=(0.5, 5.5), ylim=(0, 22000))

	ax2.bar(index - bar_width, bwa_sa_before_est, bar_width, color='navy', edgecolor='black', label='BWA')
	ax2.bar(index - bar_width, bwa_sa_after_est, bar_width, color='dodgerblue', bottom = bwa_sa_before_est, edgecolor='black', label='BWA + Estimation')
	
	ax2.bar(index, sbg_sa_before_est, bar_width, color='darkgreen', edgecolor='black', label='SBG')
	ax2.bar(index, sbg_sa_after_est, bar_width, color='seagreen', bottom = sbg_sa_before_est, edgecolor='black', label='SBG + Estimation')
	
	ax2.bar(index + bar_width, vg_sa_before_est, bar_width, color='darkorange', edgecolor='black', label='VG')
	ax2.bar(index + bar_width, vg_sa_after_est, bar_width, color='orange', bottom = vg_sa_before_est, edgecolor='black', label='VG + Estimation')
	
	# ax2.title.set_text("Run time")
	ax2.set_ylabel('Run time (sec)')
	ax2.set_xticks(np.arange(1, 6))
	ax2.set_yticks(np.arange(0, 22000, 1000))
	ax2.set_xticklabels(['S1', 'S2', 'S3', 'S4', 'S5'])
	ax2.minorticks_on()
	ax2.grid(which='major', linestyle='-', linewidth='0.05', alpha=0.3, color='black')
	ax2.grid(which='minor', linestyle='--', linewidth='0.001', alpha=0.1, color='gray')
	ax2.grid(True)

	ax2.legend(loc='upper center', fontsize='x-small')
	
	plt.tight_layout()
	plt.savefig('single_baumannii.pdf', format='pdf')
	# plt.show()
