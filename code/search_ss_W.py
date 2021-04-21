#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  4 23:39:17 2021

@author: jsheen
"""

#import matplotlib.pyplot as plt
import numpy as np
import random
import math
from operator import itemgetter
from scipy.stats import wilcoxon as wilcoxon
import pandas as pd
random.seed(1)
gen = np.random.Generator(np.random.PCG64(1))
nsim = 10000
rts = [2]
overdispersions = [0.1]
effects = [0.2, 0.4]
clusters = [100]
eits = [0.02]
nsamples = [10, 50, 100]
param_sets = []
for i in rts:
    for j in overdispersions:
        for k in clusters:
            for l in effects:
                for m in eits:
                    for n in nsamples:
                        param_sets.append([i, j, k, l , m, n])
print(param_sets)

def get_power(It_It1con_It1trt, cluster_num, nsample):
    ncluster = cluster_num
    p_vals = []
    for sim_i in range(nsim):
        chosen = random.sample(list(range(len(It_It1con_It1trt))), ncluster * 2)
        St_ItIt1_con = []
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
            sampled_St_con = sampled_t_con[0]
            sampled_It_con = sampled_t_con[2]
            sampled_It1_con = gen.multivariate_hypergeometric([St1, Et1, It1, Rt1], nsample=nsample, size=1)[0][2]
            log_It1It_con = np.log((sampled_It1_con + 1) / (sampled_It_con + 1))
            St_ItIt1_con.append([sampled_St_con, log_It1It_con])
        St_ItIt1_trt = []
        for trt_dex in range(ncluster, ncluster * 2):
            St = int(It_It1con_It1trt.iloc[chosen[trt_dex], 0])
            Et = int(It_It1con_It1trt.iloc[chosen[trt_dex], 1])
            It = int(It_It1con_It1trt.iloc[chosen[trt_dex], 2])
            Rt = int(It_It1con_It1trt.iloc[chosen[trt_dex], 3])
            St1 = int(It_It1con_It1trt.iloc[chosen[trt_dex], 16])
            Et1 = int(It_It1con_It1trt.iloc[chosen[trt_dex], 17])
            It1 = int(It_It1con_It1trt.iloc[chosen[trt_dex], 18])
            Rt1 = int(It_It1con_It1trt.iloc[chosen[trt_dex], 19])
            sampled_t_trt = gen.multivariate_hypergeometric([St, Et, It, Rt], nsample=nsample, size=1)[0]
            sampled_St_trt = sampled_t_trt[0]
            sampled_It_trt = sampled_t_trt[2]
            sampled_It1_trt = gen.multivariate_hypergeometric([St1, Et1, It1, Rt1], nsample=nsample, size=1)[0][2]
            log_It1It_trt = np.log((sampled_It1_trt + 1) / (sampled_It_trt + 1))
            St_ItIt1_trt.append([sampled_St_trt, log_It1It_trt])
        St_ItIt1_con = sorted(St_ItIt1_con, key=itemgetter(0))
        St_ItIt1_trt = sorted(St_ItIt1_trt, key=itemgetter(0))
        diffs = []
        for diff_dex in range(ncluster):
            diff = St_ItIt1_con[diff_dex][1] - St_ItIt1_trt[diff_dex][1]
            diffs.append(diff)
        p_val = wilcoxon(diffs, alternative="greater")[1]
        p_vals.append(p_val)
    power = len(np.where(np.array(p_vals) <= 0.05)[0]) / len(p_vals)
    return power
def dac(It_It1con_It1trt, min_cluster, max_cluster, nsample):
    midpoint = math.ceil((min_cluster + max_cluster) / 2)
    test_power = get_power(It_It1con_It1trt, midpoint, nsample)
    print(str(min_cluster) + ", " + str(max_cluster) + ", " + str(test_power))
    if (abs(test_power - 0.8) < 0.005):
        return midpoint
    elif ((max_cluster - min_cluster) <= 1):
        return max_cluster
    elif (min_cluster == 1 and max_cluster == 1000 and test_power < 0.4):
        return -1
    elif (min_cluster == 501 and max_cluster == 1000 and test_power < 0.6):
        return -1
    else:
        if ((test_power - 0.8) > 0):
            return dac(It_It1con_It1trt, min_cluster, midpoint, nsample)
        else:
            return dac(It_It1con_It1trt, midpoint, max_cluster, nsample)        
for param_set in param_sets:
    tgt_R0 = param_set[0]
    N_cluster = param_set[2]
    k_overdispersion = param_set[1]
    effect = param_set[3]
    expected_It_N = param_set[4]
    curr_nsample = param_set[5]
    It_It1con_It1trt = pd.read_csv("/Users/jsheen/Desktop/res/" + str(tgt_R0) + "_" + str(N_cluster) + "_" + str(k_overdispersion) + "_" + str(effect) + "_" + str(expected_It_N) + ".csv", header=None)
    It_It1con_It1trt = It_It1con_It1trt.drop([3000])
    sufficient_ncluster = dac(It_It1con_It1trt, 1, 1000, curr_nsample)
    filename = "/Users/jsheen/Desktop/res/res_" + str(tgt_R0) + "_" + str(N_cluster) + "_" + str(k_overdispersion) + "_" + str(effect) + "_" + str(expected_It_N) + "_" + str(curr_nsample) + "_W.csv"
    with open(filename, 'w') as out_f:
        out_f.write(str(sufficient_ncluster))