#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 10:50:18 2020

@author: jsheen
"""

import pandas as pd
import numpy as np
from pathlib import Path
input_folder = "/Users/jsheen/SW-CRT-outbreak/NPI_study/EoN/code_output/csvs_1_2/dists/"
output_folder = "/Users/jsheen/SW-CRT-outbreak/NPI_study/EoN/code_output/csvs_1_2/power/"

beta_vals = [0.02, 0.03, 0.04, 0.05]
ncomms = [20, 40, 60, 80]
effects = [0.2, 0.3, 0.4, 0.5, 0.6]
days = [14, 21, 28]

for beta in beta_vals:
    for day in days:
        for effect in effects:
            for ncomm in ncomms:
                print("Beta val: " + str(beta))
                print("Effect val: " + str(effect))
                print("Day val: " + str(day))
                print("ncomm val: " + str(ncomm))

                sample_sizes = list(range(1000, (ncomm * 500) + 1, 1000))
                for sample_size in sample_sizes:
                    I_folder = input_folder + "beta_" + str(beta) + "/" + "effect_" + str(effect) + "/" + "intervention_" + str(day) + "/ncomm_" + str(ncomm) + "/infect/" 
                    I_R_folder = input_folder + "beta_" + str(beta) + "/" + "effect_" + str(effect) + "/" + "intervention_" + str(day) + "/ncomm_" + str(ncomm) + "/infect_recover/"
                    I_name = I_folder + "sampled_" + str(sample_size) + ".csv"
                    I = pd.read_csv(I_name, header=None, sep=',')
                    del I[101]
                    I_R_name = I_R_folder + "sampled_" + str(sample_size) + ".csv"
                    I_R = pd.read_csv(I_R_name, header=None, sep=',')
                    del I_R[101]
                    null_I_folder = input_folder + "beta_" + str(beta) + "/" + "effect_" + str(0) + "/" + "intervention_" + str(day) + "/ncomm_" + str(ncomm) + "/infect/" 
                    null_I_R_folder = input_folder + "beta_" + str(beta) + "/" + "effect_" + str(0) + "/" + "intervention_" + str(day) + "/ncomm_" + str(ncomm) + "/infect_recover/"
                    null_I_name = null_I_folder + "sampled_" + str(sample_size) + ".csv"
                    null_I = pd.read_csv(null_I_name, header=None, sep=',')
                    del null_I[101]
                    null_I_R_name = null_I_R_folder + "sampled_" + str(sample_size) + ".csv"
                    null_I_R = pd.read_csv(null_I_R_name, header=None, sep=',')
                    del null_I_R[101]
                    
                    power_I = []
                    power_I_R = []
                    for curr_t in range(min(I.shape[0], null_I.shape[0])):
                        curr_I = np.array(list(I.iloc[curr_t]))
                        curr_null_I = np.array(list(null_I.iloc[curr_t]))
                        under_null_I = np.where(curr_I <= curr_null_I[5])[0]
                        if under_null_I.size != 0:
                            curr_power_I = under_null_I[-1] / 100
                        else:
                            curr_power_I = 0
                        power_I.append(curr_power_I)
                        curr_I_R = np.array(list(I_R.iloc[curr_t]))
                        curr_null_I_R = np.array(list(null_I_R.iloc[curr_t]))
                        under_null_I_R = np.where(curr_I_R <= curr_null_I_R[5])[0]
                        if under_null_I_R.size != 0:
                            curr_power_I_R = under_null_I_R[-1] / 100
                        else:
                            curr_power_I_R = 0
                        power_I_R.append(curr_power_I_R)
                    
                    out_I_folder = output_folder + "beta_" + str(beta) + "/" + "effect_" + str(effect) + "/" + "intervention_" + str(day) + "/ncomm_" + str(ncomm) + "/infect/sampled_" + str(sample_size + 1000) + "/"
                    Path(out_I_folder).mkdir(parents=True, exist_ok=True)
                    filename = out_I_folder + "log_I_power.csv"
                    with open(filename, 'w') as out_f:
                        for dex in range(len(power_I)):
                            out_f.write(str(power_I[dex]))
                            out_f.write("\n")
                    out_I_R_folder = output_folder + "beta_" + str(beta) + "/" + "effect_" + str(effect) + "/" + "intervention_" + str(day) + "/ncomm_" + str(ncomm) + "/infect_recover/sampled_" + str(sample_size + 1000) + "/"
                    Path(out_I_R_folder).mkdir(parents=True, exist_ok=True)
                    filename = out_I_R_folder + "log_I_R_power.csv"
                    with open(filename, 'w') as out_f:
                        for dex in range(len(power_I_R)):
                            out_f.write(str(power_I_R[dex]))
                            out_f.write("\n")
                            
                        
                        
                        
                        