#!/usr/bin/env python3.8.2


import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


if __name__ == "__main__":

	(f, (ax1,ax2)) = plt.subplots(1, 2)

	## FIRST PLOT ##

	cov = np.arange(0, 55, 5)

	bwa_ab_before_est = np.array([0, 20, 35, 51, 67, 84, 99, 113, 133, 147, 161])
	bwa_ab_after_est = np.array([0, 63, 104, 144, 182, 222, 259, 295, 335, 371, 404])
	
	sbg_ab_before_est = np.array([0, 30, 28, 38, 36, 39, 53, 59, 56, 59, 71])
	sbg_ab_after_est = np.array([0, 36, 38, 52, 55, 62, 81, 92, 94, 100, 114])
		
	vg_ab_before_est = np.array([0, 95, 187, 288, 365, 453, 553, 659, 745, 841, 929])
	vg_ab_after_est = np.array([0, 124, 233, 350, 441, 544, 733, 862, 970, 1094, 1088])

	ax1.plot(cov, bwa_ab_before_est, color='navy', linestyle='--', label='BWA Aligner')
	ax1.plot(cov, bwa_ab_after_est, color='navy', label='BWA Aligner + Estimation')
	
	ax1.plot(cov, sbg_ab_before_est, color='darkgreen', linestyle='--', label='SBG Aligner')
	ax1.plot(cov, sbg_ab_after_est, color='darkgreen', label='SBG Aligner + Estimation')
	
	ax1.plot(cov, vg_ab_before_est, color='darkorange', linestyle='--', label='VG Aligner')
	ax1.plot(cov, vg_ab_after_est, color='darkorange', label='VG Aligner + Estimation')
	
	ax1.title.set_text("A.baumannii")
	
	ax1.set_xlabel('Coverage')
	ax1.set_ylabel('Run time (sec)')
	
	ax1.set(xlim=(0, 50), ylim=(0, 1200))
	ax1.set_xticks(np.arange(0, 51, 5))
	ax1.set_yticks(np.arange(0, 1201, 100))
	ax1.minorticks_on()
	
	ax1.grid(which='major', linestyle='-', linewidth='0.05', alpha=0.3, color='black')
	ax1.grid(which='minor', linestyle='--', linewidth='0.001', alpha=0.1, color='gray')
	ax1.grid(True)

	ax1.legend(loc='upper left', fontsize='x-small')

	## SECOND PLOT ##

	bwa_sa_before_est = np.array([0, 10, 16, 24, 31, 38, 46, 52, 59, 67, 74])
	bwa_sa_after_est = np.array([0, 28, 43, 60, 76, 91, 108, 122, 137, 152, 168])
	
	sbg_sa_before_est = np.array([0, 18, 20, 22, 26, 27, 29, 30, 33, 36, 36])
	sbg_sa_after_est = np.array([0, 36, 48, 57, 69, 79, 89, 99, 109, 121, 128])
	
	vg_sa_before_est = np.array([0, 76, 128, 185, 240, 298, 352, 407, 462, 516, 624])
	vg_sa_after_est = np.array([0, 94, 156, 222, 285, 351, 414, 477, 540, 603, 720])

	ax2.plot(cov, bwa_sa_before_est, color='navy', linestyle='--', label='BWA Aligner')
	ax2.plot(cov, bwa_sa_after_est, color='navy', label='BWA Aligner + Estimation')

	ax2.plot(cov, sbg_sa_before_est, color='darkgreen', linestyle='--', label='SBG Aligner')
	ax2.plot(cov, sbg_sa_after_est, color='darkgreen', label='SBG Aligner + Estimation')
	
	ax2.plot(cov, vg_sa_before_est, color='darkorange', linestyle='--', label='VG Aligner')
	ax2.plot(cov, vg_sa_after_est, color='darkorange', label='VG Aligner + Estimation')
	
	ax2.title.set_text("S.agalactiae")
	
	ax2.set_xlabel('Coverage')
	ax2.set_ylabel('Run time (sec)')
	
	ax2.set(xlim=(0, 50), ylim=(0, 1200))
	ax2.set_xticks(np.arange(0, 51, 5))
	ax2.set_yticks(np.arange(0, 1201, 100))
	ax2.minorticks_on()
	
	ax2.grid(which='major', linestyle='-', linewidth='0.05', alpha=0.3, color='black')
	ax2.grid(which='minor', linestyle='--', linewidth='0.001', alpha=0.1, color='gray')
	ax2.grid(True)

	ax2.legend(loc='upper left', fontsize='x-small')
	
	plt.tight_layout()
	plt.savefig('coverage_experiment.pdf', format='pdf')
	# plt.show()
