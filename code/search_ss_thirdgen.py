#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 10:19:38 2021

@author: jsheen
"""

#import matplotlib.pyplot as plt
import numpy as np
import random
import math
import pandas as pd
random.seed(1)
gen = np.random.Generator(np.random.PCG64(1))
nsim = 10000
rts = [1.5, 2]
overdispersions = [0.7]
effects = [0.2]
clusters = [1000]
eits = [0.005]
nsamples = [100, 1000]
param_sets = []
for i in rts:
    for j in overdispersions:
        for k in clusters:
            for l in effects:
                for m in eits:
                    for n in nsamples:
                        param_sets.append([i, j, k, l , m, n])
def get_power(It_It1con_It1trt, cluster_num, nsample):
    ncluster = cluster_num
    log_It1It_trt = []
    log_It1It_con = []
    for sim_i in range(nsim):
        chosen = random.sample(list(range(len(It_It1con_It1trt))), ncluster * 2)
        sum_It_con = 0
        sum_It1_con = 0
        for con_dex in range(ncluster):
            St = int(It_It1con_It1trt.iloc[chosen[con_dex], 0])
            Et = int(It_It1con_It1trt.iloc[chosen[con_dex], 1])
            It = int(It_It1con_It1trt.iloc[chosen[con_dex], 2])
            Rt = int(It_It1con_It1trt.iloc[chosen[con_dex], 3])
            St1 = int(It_It1con_It1trt.iloc[chosen[con_dex], 12])
            Et1 = int(It_It1con_It1trt.iloc[chosen[con_dex], 13])
            It1 = int(It_It1con_It1trt.iloc[chosen[con_dex], 14])
            Rt1 = int(It_It1con_It1trt.iloc[chosen[con_dex], 15])
            sampled_It_con = gen.multivariate_hypergeometric([St, Et, It, Rt], nsample=nsample, size=1)[0][2]
            sampled_It1_con = gen.multivariate_hypergeometric([St1, Et1, It1, Rt1], nsample=nsample, size=1)[0][2]
            sum_It_con += sampled_It_con
            sum_It1_con += sampled_It1_con
        sum_It_trt = 0
        sum_It1_trt = 0
        for trt_dex in range(ncluster, ncluster * 2):
            St = int(It_It1con_It1trt.iloc[chosen[trt_dex], 0])
            Et = int(It_It1con_It1trt.iloc[chosen[trt_dex], 1])
            It = int(It_It1con_It1trt.iloc[chosen[trt_dex], 2])
            Rt = int(It_It1con_It1trt.iloc[chosen[trt_dex], 3])
            St1 = int(It_It1con_It1trt.iloc[chosen[trt_dex], 24])
            Et1 = int(It_It1con_It1trt.iloc[chosen[trt_dex], 25])
            It1 = int(It_It1con_It1trt.iloc[chosen[trt_dex], 26])
            Rt1 = int(It_It1con_It1trt.iloc[chosen[trt_dex], 27])
            sampled_It_trt = gen.multivariate_hypergeometric([St, Et, It, Rt], nsample=nsample, size=1)[0][2]
            sampled_It1_trt = gen.multivariate_hypergeometric([St1, Et1, It1, Rt1], nsample=nsample, size=1)[0][2]
            sum_It_trt += sampled_It_trt
            sum_It1_trt += sampled_It1_trt
        log_It1It_trt.append(np.log((sum_It1_trt + 1) / (sum_It_trt + 1)))
        log_It1It_con.append(np.log((sum_It1_con + 1) / (sum_It_con + 1)))
    log_It1It_trt.sort()
    log_It1It_con.sort()
    p_val_x = log_It1It_con[round(0.05 * len(log_It1It_con))]
    power = len(np.where(log_It1It_trt <= p_val_x)[0]) / len(log_It1It_trt)
    #plt.hist(log_It1It_trt, bins=20)
    #plt.hist(log_It1It_con, bins=20)
    #plt.xlabel("log((I_t+1 + 1) / (I_t + 1))")
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
    filename = "/Users/jsheen/Desktop/res/res_" + str(tgt_R0) + "_" + str(N_cluster) + "_" + str(k_overdispersion) + "_" + str(effect) + "_" + str(expected_It_N) + "_" + str(curr_nsample) + "_thirdgen.csv"
    with open(filename, 'w') as out_f:
        out_f.write(str(sufficient_ncluster))