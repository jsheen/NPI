#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 27 12:31:00 2021

@author: Justin Sheen

@description: script used to test stability of identified transmission rate and
              day of intervention.

"""
# Import libraries ------------------------------------------------------------
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 20})
from pathlib import Path
home = str(Path.home())
# Create parameter sets to run ------------------------------------------------
rts = [1.5, 2]
overdispersions = [0.1, 0.4, 0.7]
effects = [0.4]
clusters = [1000, 10000]
eits = [0.005]
param_sets = []
for i in rts:
    for j in overdispersions:
        for k in clusters:
            for l in effects:
                for m in eits:
                    if not (k == 10000 and j == 0.1 and i == 1.5):
                        if (k == 10000 and j == 0.1 and i == 2):
                            param_sets.append([i, j, k, l , 0.0045])
                        else:
                            param_sets.append([i, j, k, l , m])
clusters = [100]
eits = [0.02]
for i in rts:
    for j in overdispersions:
        for k in clusters:
            for l in effects:
                for m in eits:
                    param_sets.append([i, j, k, l , m])
# Test that all the transmissions are same ------------------------------------
diffs_beta = []
diffs_t = []
for param_set in param_sets:
    tgt_R0 = param_set[0]
    N_cluster = param_set[2]
    k_overdispersion = param_set[1]
    effect = param_set[3]
    expected_It_N = param_set[4]
    It_It1con_It1trt_hE = pd.read_csv(home + "/NPI/code_output/res/" + str(tgt_R0) + "_" + str(N_cluster) + "_" + str(k_overdispersion) + "_" + str(0.4) + "_" + str(expected_It_N) + ".csv", header=None)
    b_hE = It_It1con_It1trt_hE.iloc[3000, 2]
    t_hE = It_It1con_It1trt_hE.iloc[3000, 3]
    It_It1con_It1trt_lE = pd.read_csv(home + "/NPI/code_output/res/" + str(tgt_R0) + "_" + str(N_cluster) + "_" + str(k_overdispersion) + "_" + str(0.2) + "_" + str(expected_It_N) + ".csv", header=None)
    b_lE = It_It1con_It1trt_lE.iloc[3000, 2]
    t_lE = It_It1con_It1trt_lE.iloc[3000, 3]
    diffs_beta.append((b_hE - b_lE))
    diffs_t.append((t_lE - t_hE))
# Plot of distribution of diffs in beta ---------------------------------------
plt.hist(diffs_beta, color='#1f77b4', bins=10)
plt.xlabel('Difference in transmission rate, β')
plt.ylabel('Count')
plt.title('A) Algorithm 1 stability, β')
plt.tight_layout()
plt.savefig(home + '/NPI/code_output/figs/beta_stability.png', dpi=300)
max(diffs_beta)
min(diffs_beta)
np.average(np.abs(diffs_beta))
# Plot of distribution of diffs in t ------------------------------------------
plt.hist(diffs_t, color='#1f77b4', bins=10)
plt.xlabel('Difference in day of intervention, t')
plt.ylabel('Count')
plt.title('B) Algorithm 1 stability, t')
plt.tight_layout()
plt.savefig(home + '/NPI/code_output/figs/t_stability.png', dpi=300)
max(diffs_t)
min(diffs_t)
np.average(np.abs(diffs_t))