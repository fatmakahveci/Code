#!/usr/bin/env python3.7

import matplotlib.pyplot as plt
import numpy as np
import math


if __name__ == "__main__":

	exp_dict = {}

	exp_dict["vg_exp1"] = {"c1":[0.694, 0.179, 0.127], "c2":[0.669, 0.237, 0.094], "c3":[0.498, 0.454, 0.0], "c4":[0.617, 0.205, 0.178], "c5":[0.785, 0.145, 0.069]}
	exp_dict["sbg_exp1"] = {"c1":[0.694, 0.178, 0.127], "c2":[ 0.669, 0.237, 0.094], "c3":[0.499, 0.456, 0.0], "c4":[0.621, 0.207, 0.172], "c5":[0.785, 0.145, 0.069]}
	exp_dict["bwa_exp1"] = {"c1":[0.694, 0.178, 0.127], "c2":[0.67, 0.237, 0.093], "c3":[0.501, 0.453, 0.0], "c4":[0.614, 0.211, 0.175], "c5":[0.786, 0.145, 0.069]}
	
	exp_dict["vg_exp2"] = {"c1":[0.742, 0.115, 0.061], "c2":[0.758, 0.168, 0.0], "c3":[0.616, 0.332, 0.052], "c4":[0.678, 0.244, 0.078], "c5":[0.883, 0.117, 0.0]}	
	exp_dict["sbg_exp2"] = {"c1":[0.744, 0.116, 0.061], "c2":[0.76, 0.169, 0.0], "c3":[0.617, 0.333, 0.05], "c4":[0.686, 0.241, 0.073], "c5":[0.883, 0.117, 0.0]}
	exp_dict["bwa_exp2"] = {"c1":[], "c2":[], "c3":[], "c4":[], "c5":[]}
	
	exp_dict["vg_exp3"] = {"c1":[0.726, 0.148, 0.126], "c2":[0.707, 0.204, 0.089], "c3":[0.533, 0.42, 0.0], "c4":[0.589, 0.257, 0.154], "c5":[0.802, 0.133, 0.066]}
	exp_dict["sbg_exp3"] = {"c1":[0.726, 0.148, 0.126], "c2":[0.707, 0.204, 0.089], "c3":[0.537, 0.42, 0.0], "c4":[0.591, 0.258, 0.151], "c5":[0.802, 0.133, 0.066]}
	exp_dict["bwa_exp3"] = {"c1":[], "c2":[], "c3":[], "c4":[], "c5":[]}
	
	result_dict = {}
	for exp_id, exp_result_list in exp_dict.items():
		result_dict[exp_id] = []
		for cluster_id, exp_result in exp_dict[exp_id].items():
			# standard deviation 
			result = round(math.sqrt((0.77 - exp_result[0])**2+(0.15 - exp_result[1])**2+(0.08 - exp_result[2])**2),3)
			result_dict[exp_id].append(result)

	fig, axs = plt.subplots(3, 3)
	for i in range(5):
		axs[0,0].scatter(result_dict["vg_exp1"][i], result_dict["vg_exp1"][i], color='red', marker = 'o', s=50)
		axs[0,1].scatter(result_dict["sbg_exp1"][i], result_dict["sbg_exp1"][i], color='black', marker='o',s=50)
		axs[0,2].scatter(result_dict["bwa_exp1"][i], result_dict["bwa_exp1"][i], color='green', marker='o',s=50)
		
		axs[1,0].scatter(result_dict["vg_exp2"][i], result_dict["vg_exp2"][i], color='darkblue', marker='o', s=50)
		axs[1,1].scatter(result_dict["sbg_exp2"][i], result_dict["sbg_exp2"][i], color='blue', marker='o',s=50)
		axs[1,2].scatter(result_dict["bwa_exp2"][i], result_dict["bwa_exp2"][i], color='brown', marker='o', s=50)
		
		axs[2,0].scatter(result_dict["vg_exp3"][i], result_dict["vg_exp3"][i], color='pink', marker='o', s=50)
		axs[2,1].scatter(result_dict["sbg_exp3"][i], result_dict["sbg_exp3"][i], color='darkgreen', marker='o',s=50)
		axs[2,2].scatter(result_dict["bwa_exp3"][i], result_dict["bwa_exp3"][i], color='orange', marker='o', s=50)
		
	axs[0,0].title.set_text("VG (common core genes)")
	axs[0,1].title.set_text("SBG (common core genes)")
	axs[0,2].title.set_text("BWA (common core genes)")

	axs[1,0].title.set_text("VG (different core genes)")
	axs[1,1].title.set_text("SBG (different core genes)")
	axs[1,2].title.set_text("BWA (common core genes)")

	axs[2,0].title.set_text("VG (core and accessory genes)")
	axs[2,1].title.set_text("SBG (core and accessory genes)")
	axs[2,2].title.set_text("BWA (common core genes)")

	for ax in axs.flat:
		ax.set_xticks(np.arange(0, 0.5, 0.05), minor=True)
		ax.set_yticks(np.arange(0, 0.5, 0.05), minor=True)
		ax.grid(which='major', linestyle='-', linewidth='0.05', color = "black")
		ax.grid(which='minor', linestyle='--', linewidth='0.01', color = "black")
		ax.grid(True)
	
	plt.show()
