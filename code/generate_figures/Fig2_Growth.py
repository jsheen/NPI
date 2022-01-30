#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tues Oct 12 
"""
# Import libraries, set seeds and parameters ----------------------------------
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 22})
plt.rcParams['figure.figsize'] = 20, 5
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
Rt_cons_sims = []
Rt_trts_sims = []
n_cluster = 100
for i in range(200):
    chosen_con = random.sample(list(range(len(sims_con))), n_cluster)
    chosen_trt = random.sample(list(range(len(sims_trt))), n_cluster)
    I_con_full = []
    for con_dex in range(n_cluster):
        t_con_f = sims_con[chosen_con[con_dex]][0]
        I_con_f = sims_con[chosen_con[con_dex]][3]
        next_t = 0
        to_add_row = []
        for t_dex in range(len(t_con_f)):
            if t_con_f[t_dex] >= next_t:
                to_add_row.append(I_con_f[t_dex])
                next_t += 1
        I_con_full.append(to_add_row)
    I_trt_full = []
    for trt_dex in range(n_cluster):
        t_trt_f = sims_trt[chosen_trt[trt_dex]][0]
        I_trt_f = sims_trt[chosen_trt[trt_dex]][3]
        next_t = 0
        to_add_row = []
        for t_dex in range(len(t_trt_f)):
            if t_trt_f[t_dex] >= next_t:
                to_add_row.append(I_trt_f[t_dex])
                next_t += 1
        I_trt_full.append(to_add_row)
    I_con_sums = []
    I_trt_sums = []
    for t in range(200):
        I_con_sum = 0
        I_trt_sum = 0
        for comm_num in range(n_cluster):
            if len(I_con_full[comm_num]) > t:
                I_con_sum += I_con_full[comm_num][t]
            if len(I_trt_full[comm_num]) > t:
                I_trt_sum += I_trt_full[comm_num][t]
        I_con_sums.append(I_con_sum)
        I_trt_sums.append(I_trt_sum)
    Rt_cons = []
    Rt_trts = []
    for t in range(30, 189):
        if I_con_sums[t] == 0:
            Rt_con = 0
        else:
            Rt_con = I_con_sums[t + 11] / I_con_sums[t]
        if I_trt_sums[t] == 0:
            Rt_trt = 0
        else:
            Rt_trt = I_trt_sums[t + 11] / I_trt_sums[t]
        Rt_cons.append(Rt_con)
        Rt_trts.append(Rt_trt)
    Rt_cons_sims.append(Rt_cons)
    Rt_trts_sims.append(Rt_trts)
fig = plt.figure()
ax1 = fig.add_subplot(131)
for i in range(len(Rt_cons_sims)):
    ax1.plot(list(range(len(Rt_cons_sims[i]))), Rt_cons_sims[i], color='orange', linewidth=0.5, alpha=0.5)
    ax1.plot(list(range(len(Rt_trts_sims[i]))), Rt_trts_sims[i], color='grey', linewidth=0.5, alpha=0.5)
    ax1.set_ylabel('$I_{t+1} / I_t$')
    ax1.set_xlabel('Days after intervention')
    ax1.set_xlim([0, 33])
    ax1.set_ylim([0, 2])
    ax1.axhline(y=1, color='black', linewidth=0.8, linestyle='-.')
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
Rt_cons_sims = []
Rt_trts_sims = []
n_cluster = 100
for i in range(200):
    chosen_con = random.sample(list(range(len(sims_con))), n_cluster)
    chosen_trt = random.sample(list(range(len(sims_trt))), n_cluster)
    I_con_full = []
    for con_dex in range(n_cluster):
        t_con_f = sims_con[chosen_con[con_dex]][0]
        I_con_f = sims_con[chosen_con[con_dex]][3]
        next_t = 0
        to_add_row = []
        for t_dex in range(len(t_con_f)):
            if t_con_f[t_dex] >= next_t:
                to_add_row.append(I_con_f[t_dex])
                next_t += 1
        I_con_full.append(to_add_row)
    I_trt_full = []
    for trt_dex in range(n_cluster):
        t_trt_f = sims_trt[chosen_trt[trt_dex]][0]
        I_trt_f = sims_trt[chosen_trt[trt_dex]][3]
        next_t = 0
        to_add_row = []
        for t_dex in range(len(t_trt_f)):
            if t_trt_f[t_dex] >= next_t:
                to_add_row.append(I_trt_f[t_dex])
                next_t += 1
        I_trt_full.append(to_add_row)
    I_con_sums = []
    I_trt_sums = []
    for t in range(200):
        I_con_sum = 0
        I_trt_sum = 0
        for comm_num in range(n_cluster):
            if len(I_con_full[comm_num]) > t:
                I_con_sum += I_con_full[comm_num][t]
            if len(I_trt_full[comm_num]) > t:
                I_trt_sum += I_trt_full[comm_num][t]
        I_con_sums.append(I_con_sum)
        I_trt_sums.append(I_trt_sum)
    Rt_cons = []
    Rt_trts = []
    for t in range(30, 189):
        if I_con_sums[t] == 0:
            Rt_con = 0
        else:
            Rt_con = I_con_sums[t + 11] / I_con_sums[t]
        if I_trt_sums[t] == 0:
            Rt_trt = 0
        else:
            Rt_trt = I_trt_sums[t + 11] / I_trt_sums[t]
        Rt_cons.append(Rt_con)
        Rt_trts.append(Rt_trt)
    Rt_cons_sims.append(Rt_cons)
    Rt_trts_sims.append(Rt_trts)
ax2 = fig.add_subplot(132)
for i in range(len(Rt_cons_sims)):
    ax2.plot(list(range(len(Rt_cons_sims[i]))), Rt_cons_sims[i], color='orange', linewidth=0.5, alpha=0.5)
    ax2.plot(list(range(len(Rt_trts_sims[i]))), Rt_trts_sims[i], color='grey', linewidth=0.5, alpha=0.5)
    ax2.set_xlabel('Days after intervention')
    ax2.set_xlim([0, 33])
    ax2.set_ylim([0, 2])
    ax2.axhline(y=1, color='black', linewidth=0.8, linestyle='-.')
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
Rt_cons_sims = []
Rt_trts_sims = []
n_cluster = 100
for i in range(200):
    chosen_con = random.sample(list(range(len(sims_con))), n_cluster)
    chosen_trt = random.sample(list(range(len(sims_trt))), n_cluster)
    I_con_full = []
    for con_dex in range(n_cluster):
        t_con_f = sims_con[chosen_con[con_dex]][0]
        I_con_f = sims_con[chosen_con[con_dex]][3]
        next_t = 0
        to_add_row = []
        for t_dex in range(len(t_con_f)):
            if t_con_f[t_dex] >= next_t:
                to_add_row.append(I_con_f[t_dex])
                next_t += 1
        I_con_full.append(to_add_row)
    I_trt_full = []
    for trt_dex in range(n_cluster):
        t_trt_f = sims_trt[chosen_trt[trt_dex]][0]
        I_trt_f = sims_trt[chosen_trt[trt_dex]][3]
        next_t = 0
        to_add_row = []
        for t_dex in range(len(t_trt_f)):
            if t_trt_f[t_dex] >= next_t:
                to_add_row.append(I_trt_f[t_dex])
                next_t += 1
        I_trt_full.append(to_add_row)
    I_con_sums = []
    I_trt_sums = []
    for t in range(200):
        I_con_sum = 0
        I_trt_sum = 0
        for comm_num in range(n_cluster):
            if len(I_con_full[comm_num]) > t:
                I_con_sum += I_con_full[comm_num][t]
            if len(I_trt_full[comm_num]) > t:
                I_trt_sum += I_trt_full[comm_num][t]
        I_con_sums.append(I_con_sum)
        I_trt_sums.append(I_trt_sum)
    Rt_cons = []
    Rt_trts = []
    for t in range(30, 189):
        if I_con_sums[t] == 0:
            Rt_con = 0
        else:
            Rt_con = I_con_sums[t + 11] / I_con_sums[t]
        if I_trt_sums[t] == 0:
            Rt_trt = 0
        else:
            Rt_trt = I_trt_sums[t + 11] / I_trt_sums[t]
        Rt_cons.append(Rt_con)
        Rt_trts.append(Rt_trt)
    Rt_cons_sims.append(Rt_cons)
    Rt_trts_sims.append(Rt_trts)
ax3 = fig.add_subplot(133)
for i in range(len(Rt_cons_sims)):
    ax3.plot(list(range(len(Rt_cons_sims[i]))), Rt_cons_sims[i], color='orange', linewidth=0.5, alpha=0.5)
    ax3.plot(list(range(len(Rt_trts_sims[i]))), Rt_trts_sims[i], color='grey', linewidth=0.5, alpha=0.5)
    ax3.set_xlabel('Days after intervention')
    ax3.set_xlim([0, 33])
    ax3.set_ylim([0, 2])
    ax3.axhline(y=1, color='black', linewidth=0.8, linestyle='-.')
    custom_lines = [Line2D([0], [0], color='orange', lw=4),
                    Line2D([0], [0], color='grey', lw=4)]

# Plot the legend -------------------------------------------------------------
lgd = ax2.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2)

# Save figure -----------------------------------------------------------------
fig.savefig(home + '/NPI/code_output/figs/Sheen_Fig2.tiff', dpi=600, bbox_extra_artists=(lgd,), bbox_inches='tight')
fig.savefig(home + '/NPI/code_output/figs/Sheen_Fig2.eps', dpi=600, bbox_extra_artists=(lgd,), bbox_inches='tight')
fig.savefig(home + '/NPI/code_output/figs/Sheen_Fig2.png', dpi=72, bbox_extra_artists=(lgd,), bbox_inches='tight')

        
