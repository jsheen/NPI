#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 16 15:21:35 2021
"""
# Import libraries, set seeds and parameters ----------------------------------
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 22})
plt.rcParams['figure.figsize'] = 6, 16
import random
from matplotlib.lines import Line2D
import pickle
from pathlib import Path
home = str(Path.home())
random.seed(0)
palette = dict(zip(['Control', 'Treatment'], ['orange', 'grey']))
handles = [mpl.patches.Patch(color=palette[x], label=x) for x in palette.keys()]

# Plotting top panel ----------------------------------------------------------
with open(home + '/NPI/code_output/figs/Fig1_pickles/sims_con.pickle', 'rb') as handle:
    sims_con = pickle.load(handle)
with open(home + '/NPI/code_output/figs/Fig1_pickles/sims_trt.pickle', 'rb') as handle:
    sims_trt = pickle.load(handle)
to_del_con = [] # Delete any simulations that do not have any infections at time t
for i in range(len(sims_con)):
    if sims_con[i][3] is None:
        to_del_con.append(i)
sims_con = [i for j, i in enumerate(sims_con) if j not in to_del_con]
to_del_trt = []
for i in range(len(sims_trt)):
    if sims_trt[i][3] is None:
        to_del_trt.append(i)
sims_trt = [i for j, i in enumerate(sims_trt) if j not in to_del_trt]
R_cons = [] 
R_trts = []
n_cluster = 100
for i in range(200):
    chosen_con = random.sample(list(range(len(sims_con))), n_cluster)
    chosen_trt = random.sample(list(range(len(sims_trt))), n_cluster)
    R_con_full = []
    for con_dex in range(n_cluster):
        t_con_f = sims_con[chosen_con[con_dex]][0]
        R_con_f = sims_con[chosen_con[con_dex]][4]
        next_t = 0
        to_add_row = []
        for t_dex in range(len(t_con_f)):
            if t_con_f[t_dex] >= next_t:
                to_add_row.append(R_con_f[t_dex])
                next_t += 1
        R_con_full.append(to_add_row)
    R_trt_full = []
    for trt_dex in range(n_cluster):
        t_trt_f = sims_trt[chosen_trt[trt_dex]][0]
        R_trt_f = sims_trt[chosen_trt[trt_dex]][4]
        next_t = 0
        to_add_row = []
        for t_dex in range(len(t_trt_f)):
            if t_trt_f[t_dex] >= next_t:
                to_add_row.append(R_trt_f[t_dex])
                next_t += 1
        R_trt_full.append(to_add_row)
    R_con_sums = []
    R_trt_sums = []
    for t in range(200):
        R_con_sum = 0
        R_trt_sum = 0
        for comm_num in range(n_cluster):
            if len(R_con_full[comm_num]) > t:
                R_con_sum += R_con_full[comm_num][t]
            if len(R_trt_full[comm_num]) > t:
                R_trt_sum += R_trt_full[comm_num][t]
        R_con_sums.append(R_con_sum)
        R_trt_sums.append(R_trt_sum)
    R_cons.append(R_con_sums)
    R_trts.append(R_trt_sums)
fig = plt.figure()
ax1 = fig.add_subplot(311)
for i in range(len(R_cons)):
    newList = [(x / (n_cluster * 1000)) for x in R_cons[i]]
    newList2 = [(x / (n_cluster * 1000)) for x in R_trts[i]]
    ax1.plot(list(range(len(R_cons[i]))), newList, color='orange', linewidth=0.5, alpha=0.5)
    ax1.plot(list(range(len(R_trts[i]))), newList2, color='grey', linewidth=0.5, alpha=0.5)
    ax1.set_ylabel('cum. incidence')
    ax1.set_xlim([11, 63])
    ax1.set_ylim([0, 0.07])
    ax1.axvline(x=41, color='black', linewidth=0.8, linestyle='-.')
    ax1.axvline(x=30, color='black', linewidth=0.8, linestyle='-.')
    custom_lines = [Line2D([0], [0], color='orange', lw=4),
                    Line2D([0], [0], color='grey', lw=4)]
# Plotting middle panel -------------------------------------------------------
with open(home + '/NPI/code_output/figs/Fig1_pickles/sims_con_lowE.pickle', 'rb') as handle:
    sims_con = pickle.load(handle)
with open(home + '/NPI/code_output/figs/Fig1_pickles/sims_trt_lowE.pickle', 'rb') as handle:
    sims_trt = pickle.load(handle)
to_del_con = [] # Delete any simulations that do not have any infections at time t
for i in range(len(sims_con)):
    if sims_con[i][3] is None:
        to_del_con.append(i)
sims_con = [i for j, i in enumerate(sims_con) if j not in to_del_con]
to_del_trt = []
for i in range(len(sims_trt)):
    if sims_trt[i][3] is None:
        to_del_trt.append(i)
sims_trt = [i for j, i in enumerate(sims_trt) if j not in to_del_trt]
R_cons = [] 
R_trts = []
n_cluster = 100
for i in range(200):
    chosen_con = random.sample(list(range(len(sims_con))), n_cluster)
    chosen_trt = random.sample(list(range(len(sims_trt))), n_cluster)
    R_con_full = []
    for con_dex in range(n_cluster):
        t_con_f = sims_con[chosen_con[con_dex]][0]
        R_con_f = sims_con[chosen_con[con_dex]][4]
        next_t = 0
        to_add_row = []
        for t_dex in range(len(t_con_f)):
            if t_con_f[t_dex] >= next_t:
                to_add_row.append(R_con_f[t_dex])
                next_t += 1
        R_con_full.append(to_add_row)
    R_trt_full = []
    for trt_dex in range(n_cluster):
        t_trt_f = sims_trt[chosen_trt[trt_dex]][0]
        R_trt_f = sims_trt[chosen_trt[trt_dex]][4]
        next_t = 0
        to_add_row = []
        for t_dex in range(len(t_trt_f)):
            if t_trt_f[t_dex] >= next_t:
                to_add_row.append(R_trt_f[t_dex])
                next_t += 1
        R_trt_full.append(to_add_row)
    R_con_sums = []
    R_trt_sums = []
    for t in range(200):
        R_con_sum = 0
        R_trt_sum = 0
        for comm_num in range(n_cluster):
            if len(R_con_full[comm_num]) > t:
                R_con_sum += R_con_full[comm_num][t]
            if len(R_trt_full[comm_num]) > t:
                R_trt_sum += R_trt_full[comm_num][t]
        R_con_sums.append(R_con_sum)
        R_trt_sums.append(R_trt_sum)
    R_cons.append(R_con_sums)
    R_trts.append(R_trt_sums)
ax2 = fig.add_subplot(312)
for i in range(len(R_cons)):
    newList = [(x / (n_cluster * 1000)) for x in R_cons[i]]
    newList2 = [(x / (n_cluster * 1000)) for x in R_trts[i]]
    ax2.plot(list(range(len(R_cons[i]))), newList, color='orange', linewidth=0.5, alpha=0.5)
    ax2.plot(list(range(len(R_trts[i]))), newList2, color='grey', linewidth=0.5, alpha=0.5)
    ax2.set_ylabel('cum. incidence')
    ax2.set_xlim([11, 63])
    ax2.set_ylim([0, 0.07])
    ax2.axvline(x=41, color='black', linewidth=0.8, linestyle='-.')
    ax2.axvline(x=30, color='black', linewidth=0.8, linestyle='-.')
    custom_lines = [Line2D([0], [0], color='orange', lw=4),
                    Line2D([0], [0], color='grey', lw=4)]

# Plotting bottom panel -------------------------------------------------------
with open(home + '/NPI/code_output/figs/Fig1_pickles/sims_con_lowE_lowk.pickle', 'rb') as handle:
    sims_con = pickle.load(handle)
with open(home + '/NPI/code_output/figs/Fig1_pickles/sims_trt_lowE_lowk.pickle', 'rb') as handle:
    sims_trt = pickle.load(handle)
to_del_con = [] # Delete any simulations that do not have any infections at time t
for i in range(len(sims_con)):
    if sims_con[i][3] is None:
        to_del_con.append(i)
sims_con = [i for j, i in enumerate(sims_con) if j not in to_del_con]
to_del_trt = []
for i in range(len(sims_trt)):
    if sims_trt[i][3] is None:
        to_del_trt.append(i)
sims_trt = [i for j, i in enumerate(sims_trt) if j not in to_del_trt]
R_cons = [] # Create 100 simulations, randomly picking n_cluster number of clusters per arm
R_trts = []
n_cluster = 100
for i in range(200):
    chosen_con = random.sample(list(range(len(sims_con))), n_cluster)
    chosen_trt = random.sample(list(range(len(sims_trt))), n_cluster)
    R_con_full = []
    for con_dex in range(n_cluster):
        t_con_f = sims_con[chosen_con[con_dex]][0]
        R_con_f = sims_con[chosen_con[con_dex]][4]
        next_t = 0
        to_add_row = []
        for t_dex in range(len(t_con_f)):
            if t_con_f[t_dex] >= next_t:
                to_add_row.append(R_con_f[t_dex])
                next_t += 1
        R_con_full.append(to_add_row)
    R_trt_full = []
    for trt_dex in range(n_cluster):
        t_trt_f = sims_trt[chosen_trt[trt_dex]][0]
        R_trt_f = sims_trt[chosen_trt[trt_dex]][4]
        next_t = 0
        to_add_row = []
        for t_dex in range(len(t_trt_f)):
            if t_trt_f[t_dex] >= next_t:
                to_add_row.append(R_trt_f[t_dex])
                next_t += 1
        R_trt_full.append(to_add_row)
    R_con_sums = []
    R_trt_sums = []
    for t in range(200):
        R_con_sum = 0
        R_trt_sum = 0
        for comm_num in range(n_cluster):
            if len(R_con_full[comm_num]) > t:
                R_con_sum += R_con_full[comm_num][t]
            if len(R_trt_full[comm_num]) > t:
                R_trt_sum += R_trt_full[comm_num][t]
        R_con_sums.append(R_con_sum)
        R_trt_sums.append(R_trt_sum)
    R_cons.append(R_con_sums)
    R_trts.append(R_trt_sums)
ax3 = fig.add_subplot(313)
for i in range(len(R_cons)):
    newList = [(x / (n_cluster * 1000)) for x in R_cons[i]]
    newList2 = [(x / (n_cluster * 1000)) for x in R_trts[i]]
    ax3.plot(list(range(len(R_cons[i]))), newList, color='orange', linewidth=0.5, alpha=0.5)
    ax3.plot(list(range(len(R_trts[i]))), newList2, color='grey', linewidth=0.5, alpha=0.5)
    ax3.set_ylabel('cum. incidence')
    ax3.set_xlabel('Days since initial infection')
    ax3.set_xlim([11, 63])
    ax3.set_ylim([0, 0.07])
    ax3.axvline(x=41, color='black', linewidth=0.8, linestyle='-.')
    ax3.axvline(x=30, color='black', linewidth=0.8, linestyle='-.')
    custom_lines = [Line2D([0], [0], color='orange', lw=4),
                    Line2D([0], [0], color='grey', lw=4)]

# Plot the legend -------------------------------------------------------------
lgd = fig.legend(handles=handles, loc='lower center')

# Save figure -----------------------------------------------------------------
fig.savefig(home + '/NPI/code_output/figs/Fig1.tiff', dpi=600, bbox_extra_artists=(lgd,), bbox_inches='tight')
fig.savefig(home + '/NPI/code_output/figs/Fig1.png', dpi=72, bbox_extra_artists=(lgd,), bbox_inches='tight')

        
