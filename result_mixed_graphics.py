#!/usr/bin/env python3.8.2


import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


if __name__ == "__main__":

	(f, (ax1,ax2)) = plt.subplots(1, 2)

	## FIRST PLOT ##

	n_clusters = 5
	
	bwa_content1 = np.array([57.6, 64.8, 75.2, 79.1, 55.9])
	sbg_content1 = np.array([58.3, 64.6, 75.3, 79.8, 56.0])
	vg_content1 = np.array([58.1, 63.4, 75.7, 76.6, 54.4])

	bwa_content2 = np.array([23.7, 23.5, 19.3, 10.9, 23.1])
	sbg_content2 = np.array([23.2, 23.5, 19.4, 10.7, 23.3])
	vg_content2 = np.array([22.7, 24.6, 18.8, 9.5, 24.3])

	bwa_content3 = np.array([18.7, 11.7, 5.6, 9.9, 21.0])
	sbg_content3 = np.array([18.5, 11.9, 5.30, 9.60, 20.6])
	vg_content3 = np.array([19.2, 12.0, 5.5, 9.7, 21.3])

	ax1.set(xlim=(0.5, 5.5), ylim=(0, 100))

	ax1.axhline(y=66.7, linestyle='--', alpha=0.7, color='black')
	ax1.axhline(y=86.7, linestyle='--', alpha=0.7, color='black')
	ax1.axhline(y=100.0, linestyle='--', alpha=0.7, color='black')

	index = np.arange(1, n_clusters + 1)
	
	bar_width = 0.2

	ax1.set(xlim=(0.5, 5.5), ylim=(0, 100))

	ax1.bar(index - bar_width, bwa_content1, bar_width, color='navy', edgecolor='black')
	ax1.bar(index - bar_width, bwa_content2, bar_width, color='dodgerblue', bottom = bwa_content1, edgecolor='black')
	ax1.bar(index - bar_width, bwa_content3, bar_width, color='steelblue', bottom = bwa_content1+bwa_content2, edgecolor='black')

	ax1.bar(index, sbg_content1, bar_width, color='darkgreen', edgecolor='black')
	ax1.bar(index, sbg_content2, bar_width, color='seagreen', bottom = sbg_content1, edgecolor='black')
	ax1.bar(index, sbg_content3, bar_width, color='limegreen', bottom = sbg_content1+sbg_content2, edgecolor='black')

	ax1.bar(index + bar_width, vg_content1, bar_width, color='darkorange', edgecolor='black')
	ax1.bar(index + bar_width, vg_content2, bar_width, color='orange', bottom = vg_content1, edgecolor='black')
	ax1.bar(index + bar_width, vg_content3, bar_width, color='gold', bottom = vg_content1+vg_content2, edgecolor='black')

	for i in range(1,6):
		ax1.text(i-0.2, 1.0, 'BWA', color='white', fontweight='bold',ha='center', va='bottom', fontsize='x-small', rotation=90)
		ax1.text(i, 9.0, 'SBG', color='white',fontweight='bold',ha='center', va='bottom', fontsize='x-small', rotation=90)
		ax1.text(i+0.2, 16.0, 'VG', color='white',fontweight='bold',ha='center', va='bottom', fontsize='x-small', rotation=90)

	# ax1.legend(bbox_to_anchor=(1.05, 1.0), fontsize='small', loc='upper left')
	# ax1.title.set_text("Strain Content\n(Mixed Sample: 10:3:2)")
	# ax1.set_xlabel('Experiments')
	ax1.set_ylabel('Content (%)')
	ax1.set_xticks(np.arange(1, 6))
	ax1.set_yticks(np.arange(0, 110, 10))
	ax1.set_xticklabels(['E1', 'E2', 'E3', 'E4', 'E5'])
	ax1.minorticks_on()
	ax1.grid(which='major', linestyle='-', linewidth='0.05', alpha=0.3, color='black')
	ax1.grid(which='minor', linestyle='--', linewidth='0.001', alpha=0.1, color='gray')
	ax1.grid(True)

	## SECOND PLOT ##
	
	# sample = np.arange(0, 6)
	
	bwa_sa_before_est = np.array([40, 66, 57, 31, 32])
	bwa_sa_after_est = np.array([94, 127, 121, 83, 80])
	
	sbg_sa_before_est = np.array([26, 32, 32, 21, 31])
	sbg_sa_after_est = np.array([216, 284, 346, 384, 333])
		
	vg_sa_before_est = np.array([227, 239, 234, 220, 221])
	vg_sa_after_est = np.array([280, 300, 297, 271, 266])

	index = np.arange(1, n_clusters + 1)
	
	bar_width = 0.2

	ax2.set(xlim=(0.5, 5.5), ylim=(0, 800))

	ax2.bar(index - bar_width, bwa_sa_before_est, bar_width, color='navy', edgecolor='black', label='BWA')
	ax2.bar(index - bar_width, bwa_sa_after_est, bar_width, color='dodgerblue', bottom = bwa_sa_before_est, edgecolor='black', label='BWA + Estimation')
	
	ax2.bar(index, sbg_sa_before_est, bar_width, color='darkgreen', edgecolor='black', label='SBG')
	ax2.bar(index, sbg_sa_after_est, bar_width, color='seagreen', bottom = sbg_sa_before_est, edgecolor='black', label='SBG + Estimation')
	
	ax2.bar(index + bar_width, vg_sa_before_est, bar_width, color='darkorange', edgecolor='black', label='VG')
	ax2.bar(index + bar_width, vg_sa_after_est, bar_width, color='orange', bottom = vg_sa_before_est, edgecolor='black', label='VG + Estimation')
	
	# ax2.title.set_text("Run time")
	ax2.set_ylabel('Run time (sec)')
	ax2.set_xticks(np.arange(1, 6))
	ax2.set_yticks(np.arange(0, 800, 50))
	ax2.set_xticklabels(['E1', 'E2', 'E3', 'E4', 'E5'])
	ax2.minorticks_on()
	ax2.grid(which='major', linestyle='-', linewidth='0.05', alpha=0.3, color='black')
	ax2.grid(which='minor', linestyle='--', linewidth='0.001', alpha=0.1, color='gray')
	ax2.grid(True)

	ax2.legend(loc='upper right', fontsize='x-small')
	
	plt.tight_layout()
	plt.savefig('mixed_10_3_2.pdf', format='pdf')
	# plt.show()
