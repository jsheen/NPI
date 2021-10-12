#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 27 12:31:00 2021
@description: script used to test stability of identified transmission rate,
              day of intervention, and algorithm 2.
"""
# Import libraries ------------------------------------------------------------
import pandas as pd
import numpy as np
import pickle
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 25})
plt.rcParams['figure.figsize'] = 20, 5
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
    print(b_hE)
    print(b_lE)
    print((b_hE - b_lE))
    diffs_beta.append((b_hE - b_lE))
    diffs_t.append((t_lE - t_hE))
# Plot of distribution of diffs in beta ---------------------------------------
fig = plt.figure()
ax1 = fig.add_subplot(131)
ax1.hist(np.abs(diffs_beta), color='#1f77b4', bins=7)
ax1.set_xlabel('Abs. diff. in transmission rate, β')
ax1.set_ylabel('Count')
ax1.set_title('A) Algorithm 1 stability, β')
max(diffs_beta)
min(diffs_beta)
np.average(np.abs(diffs_beta))
# Plot of distribution of diffs in t ------------------------------------------
ax2 = fig.add_subplot(132)
ax2.hist(np.abs(diffs_t), color='#1f77b4', bins=7)
ax2.set_xlabel('Abs. diff. in day of intervention, t')
ax2.set_ylabel('Count')
ax2.set_title('B) Algorithm 1 stability, t')
max(diffs_t)
min(diffs_t)
np.average(np.abs(diffs_t))

# Create parameter sets to run for algorithm 2 stability ----------------------
rts = [1.5, 2]
overdispersions = [0.1, 0.4, 0.7]
effects = [0.4]
clusters = [1000, 10000]
eits = [0.005]
samplings = [100, 1000]
param_sets = []
for i in rts:
    for j in overdispersions:
        for k in clusters:
            for l in effects:
                for m in eits:
                    for n in samplings:
                        if not (k == 10000 and j == 0.1 and i == 1.5):
                            if (k == 10000 and j == 0.1 and i == 2):
                                param_sets.append([i, j, k, l , 0.0045, n])
                            else:
                                param_sets.append([i, j, k, l , m, n])
clusters = [100]
eits = [0.02]
samplings = [10, 50, 100]
for i in rts:
    for j in overdispersions:
        for k in clusters:
            for l in effects:
                for m in eits:
                    for n in samplings:
                        param_sets.append([i, j, k, l , m, n])
# Plot distribution of powers -------------------------------------------------
powers = []
for param_set in param_sets:
    tgt_R0 = param_set[0]
    N_cluster = param_set[2]
    k_overdispersion = param_set[1]
    effect = param_set[3]
    expected_It_N = param_set[4]
    sampling = param_set[5]
    # ttest -------------------------------------------------------------------
    with open(home + "/NPI/code_output/res/res_" + str(tgt_R0) + "_" + str(N_cluster) + "_" + str(k_overdispersion) + "_" + str(0.4) + "_" + str(expected_It_N) + "_" + str(sampling) + "_ttest_rerun.pickle", 'rb') as handle:
        pickle_res = pickle.load(handle)
    if len(pickle_res) > 2 and not (float(pickle_res[-1].split(',')[0]) == 999 and float(pickle_res[-1].split(',')[1]) == 1000):
        powers.append(float(pickle_res[-1].split(',')[-1]) * 100)
    # M_S ---------------------------------------------------------------------
    with open(home + "/NPI/code_output/res/res_" + str(tgt_R0) + "_" + str(N_cluster) + "_" + str(k_overdispersion) + "_" + str(0.4) + "_" + str(expected_It_N) + "_" + str(sampling) + "_M_S.pickle", 'rb') as handle:
        pickle_res = pickle.load(handle)
    if len(pickle_res) > 2 and not (float(pickle_res[-1].split(',')[0]) == 999 and float(pickle_res[-1].split(',')[1]) == 1000):
        powers.append(float(pickle_res[-1].split(',')[-1]) * 100)
    # M_SER ---------------------------------------------------------------------
    with open(home + "/NPI/code_output/res/res_" + str(tgt_R0) + "_" + str(N_cluster) + "_" + str(k_overdispersion) + "_" + str(0.4) + "_" + str(expected_It_N) + "_" + str(sampling) + "_M_SER.pickle", 'rb') as handle:
        pickle_res = pickle.load(handle)
    if len(pickle_res) > 2 and not (float(pickle_res[-1].split(',')[0]) == 999 and float(pickle_res[-1].split(',')[1]) == 1000):
        powers.append(float(pickle_res[-1].split(',')[-1]) * 100)
ax3 = fig.add_subplot(133)
ax3.hist(powers, bins=25)
ax3.set_xlabel('Power (%)')
ax3.set_ylabel('Count')
ax3.set_title('C) Algorithm 2 stability')
plt.tight_layout()
fig.savefig(home + '/NPI/code_output/figs/Sheen_FigS1.tiff', dpi=600, bbox_inches='tight')
fig.savefig(home + '/NPI/code_output/figs/Sheen_FigS1.png', dpi=72, bbox_inches='tight')
    
                        