#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 14:12:02 2021
@description: script used to create figure of the distribution of infectious individuals at time of intervention
"""
# Import libraries ------------------------------------------------------------
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 25})
plt.rcParams['figure.figsize'] = 20, 5
params = {'mathtext.default': 'regular' }          
plt.rcParams.update(params)
from pathlib import Path
home = str(Path.home())
# Create parameter sets to run ------------------------------------------------
rts = [1.5]
overdispersions = [0.7]
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
# Create figure ---------------------------------------------------------------
fig = plt.figure()

ax3 = fig.add_subplot(131)
param_set = param_sets[2]
tgt_R0 = param_set[0]
N_cluster = param_set[2]
k_overdispersion = param_set[1]
effect = param_set[3]
expected_It_N = param_set[4]
It_It1con_It1trt = pd.read_csv(home + "/NPI/code_output/res/" + str(tgt_R0) + "_" + str(N_cluster) + "_" + str(k_overdispersion) + "_" + str(effect) + "_" + str(expected_It_N) + ".csv", header=None)
It_It1con_It1trt = It_It1con_It1trt.drop([3000])
Is = It_It1con_It1trt.iloc[:, 2]
ax3.hist(Is, bins=9)
ax3.set_ylabel('Count')  
ax3.set_xlabel('$I_{t}$')
panel_letter = None
if N_cluster == 1000:
    panel_letter = 'B)'
elif N_cluster == 100:
    panel_letter = 'A)'
else:
    panel_letter = 'C)'
ax3.set_title(panel_letter + ' n=' + str(N_cluster))

ax1 = fig.add_subplot(132)
param_set = param_sets[0]
tgt_R0 = param_set[0]
N_cluster = param_set[2]
k_overdispersion = param_set[1]
effect = param_set[3]
expected_It_N = param_set[4]
It_It1con_It1trt = pd.read_csv(home + "/NPI/code_output/res/" + str(tgt_R0) + "_" + str(N_cluster) + "_" + str(k_overdispersion) + "_" + str(effect) + "_" + str(expected_It_N) + ".csv", header=None)
It_It1con_It1trt = It_It1con_It1trt.drop([3000])
Is = It_It1con_It1trt.iloc[:, 2]
ax1.hist(Is)
ax1.set_ylabel('Count')  
ax1.set_xlabel('$I_{t}$')
panel_letter = None
if N_cluster == 1000:
    panel_letter = 'B)'
elif N_cluster == 100:
    panel_letter = 'A)'
else:
    panel_letter = 'C)'
ax1.set_title(panel_letter + ' n=' + str(N_cluster))

ax2 = fig.add_subplot(133)
param_set = param_sets[1]
tgt_R0 = param_set[0]
N_cluster = param_set[2]
k_overdispersion = param_set[1]
effect = param_set[3]
expected_It_N = param_set[4]
It_It1con_It1trt = pd.read_csv(home + "/NPI/code_output/res/" + str(tgt_R0) + "_" + str(N_cluster) + "_" + str(k_overdispersion) + "_" + str(effect) + "_" + str(expected_It_N) + ".csv", header=None)
It_It1con_It1trt = It_It1con_It1trt.drop([3000])
Is = It_It1con_It1trt.iloc[:, 2]
ax2.hist(Is)
ax2.set_ylabel('Count')  
ax2.set_xlabel('$I_{t}$')
panel_letter = None
if N_cluster == 1000:
    panel_letter = 'B)'
elif N_cluster == 100:
    panel_letter = 'A)'
else:
    panel_letter = 'C)'
ax2.set_title(panel_letter + ' n=' + str(N_cluster))

plt.tight_layout()
plt.savefig(home + '/NPI/code_output/figs/FigS2.tiff', dpi=600, bbox_inches='tight')
plt.savefig(home + '/NPI/code_output/figs/FigS2.png', dpi=72, bbox_inches='tight')
