#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 17:09:02 2021

@author: jsheen
"""
import numpy as np
import pandas as pd
import os
import random
from pathlib import Path
input_folder = "/Users/jsheen/SW-CRT-outbreak/NPI_study/EoN/code_output/csvs_1_2/csvs/"
output_folder = "/Users/jsheen/SW-CRT-outbreak/NPI_study/EoN/code_output/csvs_1_2/It1_It/dists/"

beta_vals = [0.02, 0.03, 0.04, 0.05]
ncomms = [20, 40, 60, 80, 100]
effects = [0, 0.2, 0.4, 0.6]
days = [14, 21, 28]
ave_gen_time = 11
seed = 0
random.seed(seed) 
n_sims = 300

for day in days:
    for beta in beta_vals:
        for effect in effects:
            for ncomm in ncomms:
                print("Beta val: " + str(beta))
                print("Effect val: " + str(effect))
                print("Day val: " + str(day))
                print("ncomm val: " + str(ncomm))
                param_name = input_folder + "beta_" + str(beta) + "/" + "effect_" + str(effect) + "/" + "intervention_" + str(day) + "/"       
                trt_sims = []
                for nsim in range(n_sims):
                    complete_names = os.listdir(param_name + "sim_" + str(nsim))
                    treatment_names = [i for i in complete_names if "treatment" in i]
                    sub_treatment_names = random.choices(treatment_names, k=int(ncomm / 2))
                    trt_sim = None
                    for trt_comm_name in sub_treatment_names:
                        trt_comm = pd.read_csv(param_name + "sim_" + str(nsim) + "/" + trt_comm_name, header=None, sep=',')
                        del trt_comm[0]
                        del trt_comm[5]
                        if trt_sim is None:
                            trt_sim = trt_comm
                        else:
                            trt_sim = trt_comm.add(trt_sim, fill_value=0)
                    trt_sims.append(trt_sim)
                gen_It = []
                gen_It1 = []
                gen_It1_It = []
                gen_log_It1_It = []
                for curr_t in [day + ave_gen_time, day + (2 * ave_gen_time), day + (3 * ave_gen_time)]:
                    It_s = []
                    It1_s = []
                    It1_It_s = []
                    log_It1_It_s = []
                    for nsim in range(n_sims):
                        It = trt_sims[nsim].loc[day, :].values.tolist()[2]
                        It1 = trt_sims[nsim].loc[curr_t, :].values.tolist()[2]
                        It_s.append(It)
                        It1_s.append(It1)
                        It1_It_s.append((It1 + 1) / (It + 1))
                        log_It1_It_s.append(np.log((It1 + 1) / (It + 1)))
                    It_s.sort()
                    It1_s.sort()
                    It1_It_s.sort()
                    log_It1_It_s.sort()
                    gen_It.append(It_s)
                    gen_It1.append(It1_s)
                    gen_It1_It.append(It1_It_s)
                    gen_log_It1_It.append(log_It1_It_s)
                sample_output_folder = output_folder + "beta_" + str(beta) + "/effect_" + str(effect) + "/intervention_" + str(day) + "/ncomm_" + str(ncomm/2) + "/"
                Path(sample_output_folder).mkdir(parents=True, exist_ok=True)
                filename = sample_output_folder + "It.csv"
                with open(filename, 'w') as out_f:
                    for curr_gen in range(3):
                        dist = gen_It[curr_gen]
                        for dist_dex in range(len(dist)):
                            out_f.write(str(dist[dist_dex]))
                            out_f.write(",")
                        out_f.write("\n")
                filename = sample_output_folder + "It1.csv"
                with open(filename, 'w') as out_f:
                    for curr_gen in range(3):
                        dist = gen_It1[curr_gen]
                        for dist_dex in range(len(dist)):
                            out_f.write(str(dist[dist_dex]))
                            out_f.write(",")
                        out_f.write("\n")
                filename = sample_output_folder + "It1_It.csv"
                with open(filename, 'w') as out_f:
                    for curr_gen in range(3):
                        dist = gen_It1_It[curr_gen]
                        for dist_dex in range(len(dist)):
                            out_f.write(str(dist[dist_dex]))
                            out_f.write(",")
                        out_f.write("\n")
                filename = sample_output_folder + "log_It1_It.csv"
                with open(filename, 'w') as out_f:
                    for curr_gen in range(3):
                        dist = gen_log_It1_It[curr_gen]
                        for dist_dex in range(len(dist)):
                            out_f.write(str(dist[dist_dex]))
                            out_f.write(",")
                        out_f.write("\n")
                
                    