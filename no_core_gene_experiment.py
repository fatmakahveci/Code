#!/usr/bin/env python3.8.2


import matplotlib.pyplot as plt
import numpy as np
import math
import seaborn as sns

if __name__ == "__main__":

	n_clusters = 5

	(f, (ax1,ax2)) = plt.subplots(1, 2)

	## FIRST PLOT ##

	cov = [10, 25, 50, 100, 200, 300, 400, 500, 1000]

	bwa_single_ab_before_est = np.array([1074, 1090, 1121, 1170, 1214, 1251, 1249, 1259, 1309])
	bwa_single_ab_after_est = np.array([1120, 1140, 1224, 1321, 1495, 1695, 1805, 1954, 2487])
	
	sbg_single_ab_before_est = np.array([405, 382, 388, 388, 383, 382, 384, 398, 393])
	sbg_single_ab_after_est = np.array([441, 432, 489, 539, 658, 682, 942, 1099, 1556])
	
	vg_single_ab_before_est = np.array([692, 835, 999, 1361, 2150, 2974, 3982, 5347, 10776])
	vg_single_ab_after_est = np.array([786, 1013, 1177, 1661, 2512, 3536, 4686, 6088, 11998])

	ax1.plot(cov, bwa_single_ab_before_est, color='navy', linestyle='--', label='BWA Aligner')
	ax1.plot(cov, bwa_single_ab_after_est, color='navy', label='BWA Aligner + Estimation')
	
	ax1.plot(cov, sbg_single_ab_before_est, color='darkgreen', linestyle='--', label='SBG Aligner')
	ax1.plot(cov, sbg_single_ab_after_est, color='darkgreen', label='SBG Aligner + Estimation')
	
	ax1.plot(cov, vg_single_ab_before_est, color='darkorange', linestyle='--', label='VG Aligner')
	ax1.plot(cov, vg_single_ab_after_est, color='darkorange', label='VG Aligner + Estimation')
	
	ax1.title.set_text("A.baumannii")
	
	ax1.set_xlabel('Coverage')
	ax1.set_ylabel('Run time (sec)')
	
	ax1.set(xlim=(0, 1010), ylim=(0, 12010))
	ax1.set_xticks([5,40,80,200,300,400,500,600,1000])
	ax1.set_yticks(np.arange(0, 12010, 1010))
	ax1.set_xticklabels(cov, rotation=90, fontsize='x-small')
	ax1.set_yticklabels(np.arange(0, 12010, 1000),fontsize='x-small')
	# ax1.minorticks_on()
	
	ax1.grid(which='major', linestyle='-', linewidth='0.05', alpha=0.3, color='black')
	ax1.grid(which='minor', linestyle='--', linewidth='0.001', alpha=0.1, color='gray')
	ax1.grid(True)

	ax1.legend(loc='upper left', fontsize='x-small')

	## SECOND PLOT ##

	bwa_single_sa_before_est = np.array([395, 411, 426, 435, 444, 450, 453, 461, 456])
	bwa_single_sa_after_est = np.array([410, 438, 473, 482, 549, 639, 695, 840, 1431])
	
	sbg_single_sa_before_est = np.array([141, 142, 143, 143, 142, 140, 141, 139, 141])
	sbg_single_sa_after_est = np.array([145, 155, 168, 184, 223, 266, 301, 327, 583])
	
	vg_single_sa_before_est = np.array([231, 270, 294, 392, 680, 851, 1365, 1320, 2539])
	vg_single_sa_after_est = np.array([267, 351, 362, 501, 805, 1033, 1593, 1631, 3155])	

	ax2.plot(cov, bwa_single_sa_before_est, color='navy', linestyle='--', label='BWA Aligner')
	ax2.plot(cov, bwa_single_sa_after_est, color='navy', label='BWA Aligner + Estimation')

	ax2.plot(cov, sbg_single_sa_before_est, color='darkgreen', linestyle='--', label='SBG Aligner')
	ax2.plot(cov, sbg_single_sa_after_est, color='darkgreen', label='SBG Aligner + Estimation')
	
	ax2.plot(cov, vg_single_sa_before_est, color='darkorange', linestyle='--', label='VG Aligner')
	ax2.plot(cov, vg_single_sa_after_est, color='darkorange', label='VG Aligner + Estimation')
	
	ax2.title.set_text("S.agalactiae")
	
	ax2.set_xlabel('Coverage')
	ax2.set_ylabel('Run time (sec)')
	
	ax2.set(xlim=(0, 1010), ylim=(0, 12010))
	ax2.set_xticks([5,40,80,200,300,400,500,600,1000])
	ax2.set_yticks(np.arange(0, 12010, 1010))
	ax2.set_xticklabels(cov, rotation=90, fontsize='x-small')
	ax2.set_yticklabels(np.arange(0, 12010, 1000),fontsize='x-small')
	# ax2.minorticks_on()
	
	ax2.grid(which='major', linestyle='-', linewidth='0.05', alpha=0.3, color='black')
	ax2.grid(which='minor', linestyle='--', linewidth='0.001', alpha=0.1, color='gray')
	ax2.grid(True)

	ax2.legend(loc='upper left', fontsize='x-small')
	
	plt.tight_layout()
	plt.savefig('no_core_gene_experiment.pdf', format='pdf')
	# plt.show()
