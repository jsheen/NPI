#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 27 12:31:00 2021
@description: script used to test stability of identified transmission rate,
              day of intervention, and algorithm 2.
"""
# Import libraries ------------------------------------------------------------
import pickle
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 22})
plt.rcParams['figure.figsize'] = 7, 5
from pathlib import Path
home = str(Path.home())
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
plt.hist(powers, bins=30)
plt.xlabel('Power (%)')
plt.ylabel('Count')
plt.title('Algorithm 2 stability')
plt.tight_layout()
plt.savefig(home + '/NPI/code_output/figs/Sheen_FigS1.tiff', dpi=600, bbox_inches='tight')
plt.savefig(home + '/NPI/code_output/figs/Sheen_FigS1.png', dpi=72, bbox_inches='tight')
    
                        