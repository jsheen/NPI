#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 21:48:17 2020

@author: jsheen
"""

# Import libraries and set input and output folders ---------------------------
import pandas as pd
import numpy as np
from pathlib import Path
import os
import random
input_folder = "/Users/jsheen/SW-CRT-outbreak/NPI_study/EoN/code_output/csvs_1_2/csvs/"
output_folder = "/Users/jsheen/SW-CRT-outbreak/NPI_study/EoN/code_output/csvs_1_2/dists/"

beta_vals = [0.02, 0.03, 0.04, 0.05]
ncomms = [40, 60, 80, 100]
effects = [0, 0.2, 0.4, 0.6]
days = [14, 21, 28]

seed = 0
random.seed(seed) 
gen = np.random.Generator(np.random.PCG64(seed))

n_sims = 300

for beta in beta_vals:
    for day in days:
        for effect in effects:
            for ncomm in ncomms:
                print("Beta val: " + str(beta))
                print("Effect val: " + str(effect))
                print("Day val: " + str(day))
                print("ncomm val: " + str(ncomm))
                param_name = input_folder + "beta_" + str(beta) + "/" + "effect_" + str(effect) + "/" + "intervention_" + str(day) + "/"
                
                con_sims = []
                trt_sims = []
                for nsim in range(n_sims):
                    complete_names = os.listdir(param_name + "sim_" + str(nsim))
                    control_names = [i for i in complete_names if "control" in i]
                    treatment_names = [i for i in complete_names if "treatment" in i]
                    sub_control_names = random.choices(control_names, k=int(ncomm / 2))
                    sub_treatment_names = random.choices(treatment_names, k=int(ncomm / 2))
                    con_sim = None
                    for con_comm_name in sub_control_names:
                        con_comm = pd.read_csv(param_name + "sim_" + str(nsim) + "/" + con_comm_name, header=None, sep=',')
                        del con_comm[0]
                        del con_comm[5]
                        if con_sim is None:
                            con_sim = con_comm
                        else:
                            con_sim = con_comm.add(con_sim, fill_value=0)
                    con_sims.append(con_sim)
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
                
                min_con_sims = min([i.shape[0] for i in con_sims]) - 1
                min_trt_sims = min([i.shape[0] for i in trt_sims]) - 1
                min_con_trt_sims = min(min_con_sims, min_trt_sims)
                
                sample_sizes = list(range(1000, (int(ncomm * 500) + 1000), 1000))
                sample_sizes_I_res = []
                sample_sizes_I_R_res = []
                for curr_t in range(min_con_trt_sims + 1):
                    print("Curr. time: " + str(curr_t))
                    t_log_ratio_I_pct = []
                    t_log_ratio_I_R_pct = []
                    for sample_size in sample_sizes:
                        log_ratio_I = []
                        log_ratio_I_R = []
                        for nsim in range(n_sims):
                            curr_trt = trt_sims[nsim].loc[curr_t, :].values.tolist()
                            sampled_trt = gen.multivariate_hypergeometric(curr_trt, nsample=int(sample_size / 2), size=500)
                            curr_con = con_sims[nsim].loc[curr_t, :].values.tolist()
                            sampled_con = gen.multivariate_hypergeometric(curr_con, nsample=int(sample_size / 2), size=500)
                            log_ratio_I_sim = np.log((sampled_trt[:,2] + 1) / (sampled_con[:,2] + 1))
                            log_ratio_I_R_sim = np.log((sampled_trt[:,2] + sampled_trt[:,3] + 1) / (sampled_con[:,2] + sampled_con[:,3] + 1))
                            log_ratio_I.extend(list(log_ratio_I_sim))
                            log_ratio_I_R.extend(list(log_ratio_I_R_sim))
                        log_ratio_I.sort()
                        log_ratio_I_R.sort()
                        indices_pct = list(range(-1, len(log_ratio_I), 1500))
                        indices_pct[0] = 0
                        t_log_ratio_I_pct.append([log_ratio_I[x] for x in indices_pct])
                        t_log_ratio_I_R_pct.append([log_ratio_I_R[x] for x in indices_pct])
                    sample_sizes_I_res.append(t_log_ratio_I_pct)
                    sample_sizes_I_R_res.append(t_log_ratio_I_R_pct)
                
                for sample_size_i in range(len(sample_sizes)):
                    sample_output_folder = output_folder + "beta_" + str(beta) + "/effect_" + str(effect) + "/intervention_" + str(day) + "/ncomm_" + str(ncomm) + "/infect/"
                    Path(sample_output_folder).mkdir(parents=True, exist_ok=True)
                    filename = sample_output_folder + "sampled_" + str(int((sample_size_i + 1) * 1000)) + ".csv"
                    with open(filename, 'w') as out_f:
                        for curr_t in range(len(sample_sizes_I_res)):
                            dist_pct = sample_sizes_I_res[curr_t][sample_size_i]
                            for dist_dex in range(len(dist_pct)):
                                out_f.write(str(dist_pct[dist_dex]))
                                out_f.write(",")
                            out_f.write("\n")
                    sample_output_folder = output_folder + "beta_" + str(beta) + "/effect_" + str(effect) + "/intervention_" + str(day) + "/ncomm_" + str(ncomm) + "/infect_recover/"
                    Path(sample_output_folder).mkdir(parents=True, exist_ok=True)
                    filename = sample_output_folder + "sampled_" + str(int((sample_size_i + 1) * 1000)) + ".csv"
                    with open(filename, 'w') as out_f:
                        for curr_t in range(len(sample_sizes_I_R_res)):
                            dist_pct = sample_sizes_I_R_res[curr_t][sample_size_i]
                            for dist_dex in range(len(dist_pct)):
                                out_f.write(str(dist_pct[dist_dex]))
                                out_f.write(",")
                            out_f.write("\n")
                
                traj_folder_name = output_folder + "beta_" + str(beta) + "/effect_" + str(effect) + "/intervention_" + str(day) + "/ncomm_" + str(ncomm) + "/traj/"
                Path(traj_folder_name).mkdir(parents=True, exist_ok=True)
                filename = traj_folder_name + "control_S.csv"
                with open(filename, 'w') as out_f:
                    for con_dex in range(len(con_sims)):
                        for curr_t in range(len(con_sims[con_dex][1])):
                            out_f.write(str(con_sims[con_dex][1][curr_t]))
                            out_f.write(",")
                        out_f.write("\n")
                filename = traj_folder_name + "control_E.csv"
                with open(filename, 'w') as out_f:
                    for con_dex in range(len(con_sims)):
                        for curr_t in range(len(con_sims[con_dex][2])):
                            out_f.write(str(con_sims[con_dex][2][curr_t]))
                            out_f.write(",")
                        out_f.write("\n")
                filename = traj_folder_name + "control_I.csv"
                with open(filename, 'w') as out_f:
                    for con_dex in range(len(con_sims)):
                        for curr_t in range(len(con_sims[con_dex][3])):
                            out_f.write(str(con_sims[con_dex][3][curr_t]))
                            out_f.write(",")
                        out_f.write("\n")
                filename = traj_folder_name + "control_R.csv"
                with open(filename, 'w') as out_f:
                    for con_dex in range(len(con_sims)):
                        for curr_t in range(len(con_sims[con_dex][4])):
                            out_f.write(str(con_sims[con_dex][4][curr_t]))
                            out_f.write(",")
                        out_f.write("\n")
                filename = traj_folder_name + "treatment_S.csv"
                with open(filename, 'w') as out_f:
                    for trt_dex in range(len(trt_sims)):
                        for curr_t in range(len(trt_sims[trt_dex][1])):
                            out_f.write(str(trt_sims[trt_dex][1][curr_t]))
                            out_f.write(",")
                        out_f.write("\n")
                filename = traj_folder_name + "treatment_E.csv"
                with open(filename, 'w') as out_f:
                    for trt_dex in range(len(trt_sims)):
                        for curr_t in range(len(trt_sims[trt_dex][2])):
                            out_f.write(str(trt_sims[trt_dex][2][curr_t]))
                            out_f.write(",")
                        out_f.write("\n")
                filename = traj_folder_name + "treatment_I.csv"
                with open(filename, 'w') as out_f:
                    for trt_dex in range(len(trt_sims)):
                        for curr_t in range(len(trt_sims[trt_dex][3])):
                            out_f.write(str(trt_sims[trt_dex][3][curr_t]))
                            out_f.write(",")
                        out_f.write("\n")
                filename = traj_folder_name + "treatment_R.csv"
                with open(filename, 'w') as out_f:
                    for trt_dex in range(len(trt_sims)):
                        for curr_t in range(len(trt_sims[trt_dex][4])):
                            out_f.write(str(trt_sims[trt_dex][4][curr_t]))
                            out_f.write(",")
                        out_f.write("\n")
                
                    
                    
                
                
            
            