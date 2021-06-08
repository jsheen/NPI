#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  4 22:23:49 2021

@author: Justin
"""

# Import libraries and set seeds ----------------------------------------------
import numpy as np
import random
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
home = str(Path.home())
random.seed(1)
gen = np.random.Generator(np.random.PCG64(1))
# Create parameter sets to run ------------------------------------------------
rts = [1.5, 2]
overdispersions = [0.1, 0.4, 0.7]
clusters = [1000, 10000]
effects = [0.2, 0.4]
eits = [0.005]
param_sets = []
for i in rts:
    for j in eits:
        for k in clusters:
            for l in effects:
                for m in overdispersions:
                    if not (k == 10000 and m == 0.1 and i == 1.5):
                        if (k == 10000 and m == 0.1):
                    	    param_sets.append([i, 0.0045, k, l , m])
                        else:
                            param_sets.append([i, j, k, l , m])
clusters = [100]
eits = [0.02]
for i in rts:
    for j in eits:
        for k in clusters:
            for l in effects:
                for m in overdispersions:
                    param_sets.append([i, j, k, l , m])
# Run for each param set ------------------------------------------------------
to_pd = []
for param_set in param_sets:
    tgt_R0 = param_set[0]
    N_cluster = param_set[2]
    k_overdispersion = param_set[4]
    effect = param_set[3]
    expected_It_N = param_set[1]
    It_It1con_It1trt = pd.read_csv(home + "/NPI/code_output/res/" + str(tgt_R0) + "_" + str(N_cluster) + "_" + str(k_overdispersion) + "_" + str(effect) + "_" + str(expected_It_N) + ".csv", header=None)
    beta = It_It1con_It1trt[2][3000]
    It_It1con_It1trt = It_It1con_It1trt.drop([3000])
    beta_ests_0s = []
    beta_ests = []
    beta_est_cs = []
    cnt_0s = 0
    cnt_1s = 0
    Sts = []
    Its = []
    It1s = []
    for sim_num in range(len(It_It1con_It1trt)):
        St = It_It1con_It1trt[0][sim_num]
        It = It_It1con_It1trt[2][sim_num]
        It1 = It_It1con_It1trt[6][sim_num]
        beta_est = It1 / (St * It)
        beta_ests_0s.append(beta_est)
        if beta_est != 0:
            beta_ests.append(beta_est)
            Sts.append(St)
            Its.append(It)
            It1s.append(It1)
        else:
            cnt_0s += 1
        if It1 == 1:
            cnt_1s += 1
        beta_est_c = (It1 + 1) / (It + 1)
        beta_est_cs.append(beta_est_c)
    plt.hist(beta_ests)
    plt.axvline(beta, color='red')
    plt.title("n=" + str(N_cluster) + " k=" + str(k_overdispersion) +
              " 0s=" + str(cnt_0s) + " 1s=" + str(cnt_1s))
    plt.xlabel('estimated betas')
    plt.ylabel('count')
    plt.show()
    plt.clf()
    to_pd.append([tgt_R0, N_cluster, k_overdispersion, beta, np.mean(beta_ests), 
                  np.mean(beta_ests_0s), beta - np.mean(beta_ests), beta - np.mean(beta_ests_0s), 
                  np.mean(Sts), np.mean(Its), np.mean(It1s), np.mean(It1s) / np.mean(Its),
                  1 - (np.mean(beta_ests) / beta)])
df = pd.DataFrame(to_pd, columns=['R0', 'N', 'k', 'true_beta', 'beta_est', 'beta_est wih 0s', 
                                  'diff beta_est', 'diff beta_est with 0s', 
                                  'St', 'It', 'It1', 'R00', 'missing'])
df.to_csv(home + '/Desktop/beta_sanity.csv', index=False)
plt.scatter(df['true_beta'], df['beta_est'])        
plt.xlabel('true beta')
plt.ylabel('estimated beta')
plt.xlim(0, 0.032)
plt.ylim(0, 0.032)        
plt.plot([0,1], [0,1], 'k-', color = 'r', linestyle='dashed')
plt.show()
plt.clf()
plt.scatter(df['true_beta'], df['beta_est wih 0s'])        
plt.xlabel('true beta')
plt.ylabel('estimated beta including 0s')
plt.xlim(0, 0.032)
plt.ylim(0, 0.032)        
plt.plot([0,1], [0,1], 'k-', color = 'r', linestyle='dashed')
plt.show()
plt.clf()


                    