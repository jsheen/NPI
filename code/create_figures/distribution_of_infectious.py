#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 14:12:02 2021

@author: Justin Sheen

@description: script used to create figure of the distribution of infectious individuals at time of intervention

"""
# Import libraries ------------------------------------------------------------
import pandas as pd
import matplotlib.pyplot as plt
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
# Create parameter sets to run ------------------------------------------------
for param_set in param_sets:
    tgt_R0 = param_set[0]
    N_cluster = param_set[2]
    k_overdispersion = param_set[1]
    effect = param_set[3]
    expected_It_N = param_set[4]
    It_It1con_It1trt = pd.read_csv(home + "/NPI/res/" + str(tgt_R0) + "_" + str(N_cluster) + "_" + str(k_overdispersion) + "_" + str(effect) + "_" + str(expected_It_N) + ".csv", header=None)
    It_It1con_It1trt = It_It1con_It1trt.drop([3000])
    Is = It_It1con_It1trt.iloc[:, 2]
    # Plot distribution of infectious individuals at time of intervention -----
    plt.hist(Is)
    plt.ylabel('Frequency')  
    params = {'mathtext.default': 'regular' }          
    plt.rcParams.update(params)
    plt.xlabel('$I_{t}$')
    panel_letter = None
    if N_cluster == 1000:
        panel_letter = 'B)'
    elif N_cluster == 100:
        panel_letter = 'A)'
    else:
        panel_letter = 'C)'
    plt.title(panel_letter + ' n=' + str(N_cluster))
    plt.show()

