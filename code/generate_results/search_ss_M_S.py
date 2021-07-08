#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  5 11:22:06 2021

@description: algorithm to find sufficient number of clusters when matching on
              number of susceptible individuals and performing a paired one-sample 
              t-test on the differences of the outcome.

"""
# Import libraries and set seeds ----------------------------------------------
import numpy as np
import random
import math
from scipy.stats import ttest_1samp as ttest_1samp
import pandas as pd
import pickle
from pathlib import Path
home = str(Path.home())
random.seed(1)
gen = np.random.Generator(np.random.PCG64(1))
# Create parameter sets to run ------------------------------------------------
nsim = 10000
rts = [1.5, 2]
overdispersions = [0.1, 0.4, 0.7]
effects = [0.2, 0.4]
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
# Functions for algorithm ----------------------------------------------------
def get_power(It_It1con_It1trt, cluster_num, nsample):
    ncluster = cluster_num
    p_vals = []
    for sim_i in range(nsim):
        chosen = random.sample(list(range(len(It_It1con_It1trt))), ncluster * 2)
        St_ItIt1 = []
        for clus_dex in range(ncluster * 2):
            St = int(It_It1con_It1trt.iloc[chosen[clus_dex], 0])
            Et = int(It_It1con_It1trt.iloc[chosen[clus_dex], 1])
            It = int(It_It1con_It1trt.iloc[chosen[clus_dex], 2])
            Rt = int(It_It1con_It1trt.iloc[chosen[clus_dex], 3])
            St1_con = int(It_It1con_It1trt.iloc[chosen[clus_dex], 4])
            Et1_con = int(It_It1con_It1trt.iloc[chosen[clus_dex], 5])
            It1_con = int(It_It1con_It1trt.iloc[chosen[clus_dex], 6])
            Rt1_con = int(It_It1con_It1trt.iloc[chosen[clus_dex], 7])
            St1_trt = int(It_It1con_It1trt.iloc[chosen[clus_dex], 16])
            Et1_trt = int(It_It1con_It1trt.iloc[chosen[clus_dex], 17])
            It1_trt = int(It_It1con_It1trt.iloc[chosen[clus_dex], 18])
            Rt1_trt = int(It_It1con_It1trt.iloc[chosen[clus_dex], 19])
            sampled_t = gen.multivariate_hypergeometric([St, Et, It, Rt], nsample=nsample, size=1)[0]
            sampled_St = sampled_t[0]
            sampled_It = sampled_t[2]
            sampled_It1_con = gen.multivariate_hypergeometric([St1_con, Et1_con, It1_con, Rt1_con], nsample=nsample, size=1)[0][2]
            log_It1It_con = np.log((sampled_It1_con + 1) / (sampled_It + 1))
            sampled_It1_trt = gen.multivariate_hypergeometric([St1_trt, Et1_trt, It1_trt, Rt1_trt], nsample=nsample, size=1)[0][2]
            log_It1It_trt = np.log((sampled_It1_trt + 1) / (sampled_It + 1))
            St_ItIt1.append([sampled_St, log_It1It_con, log_It1It_trt])
        St_ItIt1_con_trt = []
        while len(St_ItIt1) != 0:
            dists = [abs(St_ItIt1[0][0] - x[0]) for x in St_ItIt1[1:]]
            index_clus2 = dists.index(min(dists)) + 1
            St_ItIt1_con_trt.append([St_ItIt1[0], St_ItIt1[index_clus2]])
            del St_ItIt1[index_clus2]
            del St_ItIt1[0]
        diffs = [x[0][1] - x[1][2] for x in St_ItIt1_con_trt]
        if len(diffs) != ncluster:
            raise NameError("Incorrect number of pairs.")
        p_val = ttest_1samp(diffs, popmean=0, alternative="greater")[1]
        p_vals.append(p_val)
    if len(p_vals) != nsim:
        raise NameError("Incorrect number of p_vals.")
    power = len(np.where(np.array(p_vals) <= 0.05)[0]) / len(p_vals)
    return power
def dac(It_It1con_It1trt, min_cluster, max_cluster, nsample, iters):
    midpoint = math.ceil((min_cluster + max_cluster) / 2)
    test_power = get_power(It_It1con_It1trt, midpoint, nsample)
    print(str(min_cluster) + ", " + str(max_cluster) + ", " + str(test_power))
    curr_res = str(min_cluster) + ", " + str(max_cluster) + ", " + str(test_power)
    iters.append(curr_res)
    if (abs(test_power - 0.8) < 0.005):
        return midpoint, iters
    elif ((max_cluster - min_cluster) <= 1):
        return max_cluster, iters
    elif (min_cluster == 1 and max_cluster == 1000 and test_power < 0.4):
        return -1, iters
    elif (min_cluster == 501 and max_cluster == 1000 and test_power < 0.6):
        return -1, iters
    else:
        if ((test_power - 0.8) > 0):
            return dac(It_It1con_It1trt, min_cluster, midpoint, nsample, iters)
        else:
            return dac(It_It1con_It1trt, midpoint, max_cluster, nsample, iters)
# Solve for all parameter sets ------------------------------------------------       
for param_set in param_sets:
    tgt_R0 = param_set[0]
    N_cluster = param_set[2]
    k_overdispersion = param_set[1]
    effect = param_set[3]
    expected_It_N = param_set[4]
    curr_nsample = param_set[5]
    It_It1con_It1trt = pd.read_csv(home + "/NPI/code_output/res/" + str(tgt_R0) + "_" + str(N_cluster) + "_" + str(k_overdispersion) + "_" + str(effect) + "_" + str(expected_It_N) + ".csv", header=None)
    It_It1con_It1trt = It_It1con_It1trt.drop([3000])
    sufficient_ncluster, focal_iters = dac(It_It1con_It1trt, 1, 1000, curr_nsample, [])
    filename = home + "/NPI/code_output/res/res_" + str(tgt_R0) + "_" + str(N_cluster) + "_" + str(k_overdispersion) + "_" + str(effect) + "_" + str(expected_It_N) + "_" + str(curr_nsample) + "_M_S.csv"
    with open(filename, 'w') as out_f:
        out_f.write(str(sufficient_ncluster))
    with open(home + '/NPI/code_output/res/res_' + str(tgt_R0) + "_" + str(N_cluster) + "_" + str(k_overdispersion) + "_" + str(effect) + "_" + str(expected_It_N) + "_" + str(curr_nsample) + '_M_S.pickle', 'wb') as handle:
        pickle.dump(focal_iters, handle, protocol=pickle.HIGHEST_PROTOCOL)
        
        
        
        
        