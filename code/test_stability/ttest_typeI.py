#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 10 11:41:42 2021

@author: jsheen
"""

# Import libraries and set seeds ----------------------------------------------
import numpy as np
import random
from scipy.stats import ttest_ind as ttest_ind
import pandas as pd
from pathlib import Path
home = str(Path.home())
random.seed(1)
gen = np.random.Generator(np.random.PCG64(1))
'''
# Create parameter sets to run ------------------------------------------------
nsim = 1000
rts = [2]
overdispersions = [0.7]
effects = [0.4]
clusters = [1000, 10000]
eits = [0.005]
nsamples = [100, 1000]
param_sets = []
for i in rts:
    for j in overdispersions:
        for k in clusters:
            for l in effects:
                for m in eits:
                    for n in nsamples:
                        if not (k == 10000 and j == 0.1 and i == 1.5):
                            if k == 10000 and j == 0.1 and i == 2:
                                param_sets.append([i, j, k, l , 0.0045, n])
                            else:
                                param_sets.append([i, j, k, l , m, n])
clusters = [100]
eits = [0.02]
nsamples = [10, 50, 100]
for i in rts:
    for j in overdispersions:
        for k in clusters:
            for l in effects:
                for m in eits:
                    for n in nsamples:
                        param_sets.append([i, j, k, l , m, n])
print(param_sets)
def get_error(It_It1con_It1trt, cluster_num, nsample):
    ncluster = cluster_num
    p_vals = []
    for sim_i in range(nsim):
        chosen = random.sample(list(range(len(It_It1con_It1trt))), ncluster * 2)
        ItIt1_con = []
        for con_dex in range(ncluster):
            St = int(It_It1con_It1trt.iloc[chosen[con_dex], 0])
            Et = int(It_It1con_It1trt.iloc[chosen[con_dex], 1])
            It = int(It_It1con_It1trt.iloc[chosen[con_dex], 2])
            Rt = int(It_It1con_It1trt.iloc[chosen[con_dex], 3])
            St1 = int(It_It1con_It1trt.iloc[chosen[con_dex], 4])
            Et1 = int(It_It1con_It1trt.iloc[chosen[con_dex], 5])
            It1 = int(It_It1con_It1trt.iloc[chosen[con_dex], 6])
            Rt1 = int(It_It1con_It1trt.iloc[chosen[con_dex], 7])
            sampled_t_con = gen.multivariate_hypergeometric([St, Et, It, Rt], nsample=nsample, size=1)[0]
            sampled_It_con = sampled_t_con[2]
            sampled_It1_con = gen.multivariate_hypergeometric([St1, Et1, It1, Rt1], nsample=nsample, size=1)[0][2]
            log_It1It_con = np.log((sampled_It1_con + 1) / (sampled_It_con + 1))
            ItIt1_con.append(log_It1It_con)
        ItIt1_trt = []
        for trt_dex in range(ncluster, ncluster * 2):
            St = int(It_It1con_It1trt.iloc[chosen[trt_dex], 0])
            Et = int(It_It1con_It1trt.iloc[chosen[trt_dex], 1])
            It = int(It_It1con_It1trt.iloc[chosen[trt_dex], 2])
            Rt = int(It_It1con_It1trt.iloc[chosen[trt_dex], 3])
            St1 = int(It_It1con_It1trt.iloc[chosen[trt_dex], 4])
            Et1 = int(It_It1con_It1trt.iloc[chosen[trt_dex], 5])
            It1 = int(It_It1con_It1trt.iloc[chosen[trt_dex], 6])
            Rt1 = int(It_It1con_It1trt.iloc[chosen[trt_dex], 7])
            sampled_t_trt = gen.multivariate_hypergeometric([St, Et, It, Rt], nsample=nsample, size=1)[0]
            sampled_It_trt = sampled_t_trt[2]
            sampled_It1_trt = gen.multivariate_hypergeometric([St1, Et1, It1, Rt1], nsample=nsample, size=1)[0][2]
            log_It1It_trt = np.log((sampled_It1_trt + 1) / (sampled_It_trt + 1))
            ItIt1_trt.append(log_It1It_trt)
        p_val = ttest_ind(ItIt1_con, ItIt1_trt, alternative='greater', equal_var=False)[1]
        p_vals.append(p_val)
    power = len(np.where(np.array(p_vals) <= 0.05)[0]) / len(p_vals)
    return power
for param_set in param_sets:
    print(param_set)
    tgt_R0 = param_set[0]
    N_cluster = param_set[2]
    k_overdispersion = param_set[1]
    effect = param_set[3]
    expected_It_N = param_set[4]
    curr_nsample = param_set[5]
    It_It1con_It1trt = pd.read_csv(home + "/NPI/code_output/res/" + str(tgt_R0) + "_" + str(N_cluster) + "_" + str(k_overdispersion) + "_" + str(effect) + "_" + str(expected_It_N) + ".csv", header=None)
    It_It1con_It1trt = It_It1con_It1trt.drop([3000])
    ncluster_res = pd.read_csv(home + "/NPI/code_output/res/res_" + str(tgt_R0) + "_" + str(N_cluster) + "_" + str(k_overdispersion) + "_" + str(effect) + "_" + str(expected_It_N) + "_" + str(curr_nsample) + "_ttest_rerun.csv", header=None)
    ncluster = ncluster_res[0][0]
    if ncluster != -1 and ncluster != 1000:
        to_save = get_error(It_It1con_It1trt, ncluster, curr_nsample)
        with open(home + '/Desktop/typeI/res_' + str(tgt_R0) + "_" + str(N_cluster) + "_" + str(k_overdispersion) + "_" + str(effect) + "_" + str(expected_It_N) + "_" + str(curr_nsample) + '_ttest_typeI.csv', 'w') as out_f:
            out_f.write(str(to_save))
    else:
        to_save = -1
        with open(home + '/Desktop/typeI/res_' + str(tgt_R0) + "_" + str(N_cluster) + "_" + str(k_overdispersion) + "_" + str(effect) + "_" + str(expected_It_N) + "_" + str(curr_nsample) + '_ttest_typeI.csv', 'w') as out_f:
            out_f.write(str(to_save))
'''
# Create csv of results -------------------------------------------------------
rts = [1.5, 2]
overdispersions = [0.1, 0.4, 0.7]
effects = [0.2, 0.4]
clusters = [1000, 10000]
eits = [0.005]
nsamples = [100, 1000]
param_sets = []
for i in rts:
    for j in clusters:
        for k in overdispersions:
            for l in effects:
                for m in eits:
                    for n in nsamples:
                        if not (j == 10000 and k == 0.1 and i == 1.5):
                            if j == 10000 and k == 0.1 and i == 2:
                                param_sets.append([i, j, k, l , 0.0045, n])
                            else:
                                param_sets.append([i, j, k, l , m, n])
clusters = [100]
eits = [0.02]
nsamples = [10, 50, 100]
for i in rts:
    for j in clusters:
        for k in overdispersions:
            for l in effects:
                for m in eits:
                    for n in nsamples:
                        param_sets.append([i, j, k, l , m, n])
print(param_sets)
to_pd = []
for param_set in param_sets:
    tgt_R0 = param_set[0]
    N_cluster = param_set[1]
    k_overdispersion = param_set[2]
    effect = param_set[3]
    expected_It_N = param_set[4]
    curr_nsample = param_set[5]
    param_res = pd.read_csv(home + '/Desktop/typeI/res_' + str(tgt_R0) + "_" + str(N_cluster) + "_" + str(k_overdispersion) + "_" + str(effect) + "_" + str(expected_It_N) + "_" + str(curr_nsample) + '_ttest_typeI.csv')
    new_row = []
    new_row.append(tgt_R0)
    new_row.append(N_cluster)
    new_row.append(k_overdispersion)
    new_row.append(effect)
    new_row.append(curr_nsample)
    new_row.append(param_res.columns[0])
    to_pd.append(new_row)
df = pd.DataFrame(to_pd, columns = ['R0', 'N', 'k', 'true effect', 'n sample', 'error']) 
df.to_csv(home + '/Desktop/typeI/error.csv', index=False)










