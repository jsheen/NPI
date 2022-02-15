#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 16 15:21:35 2021
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
props = dict(boxstyle='round', facecolor='white', alpha=0.5)

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
EIR_cons = [] 
EIR_trts = []
n_cluster = 100
for i in range(200):
    chosen_con = random.sample(list(range(len(sims_con))), n_cluster)
    chosen_trt = random.sample(list(range(len(sims_trt))), n_cluster)
    # Control clusters
    E_con_full = []
    I_con_full = []
    R_con_full = []
    for con_dex in range(n_cluster):
        t_con_f = sims_con[chosen_con[con_dex]][0]
        E_con_f = sims_con[chosen_con[con_dex]][2]
        I_con_f = sims_con[chosen_con[con_dex]][3]
        R_con_f = sims_con[chosen_con[con_dex]][4]
        next_t = 0
        to_add_row_con_E = []
        to_add_row_con_I = []
        to_add_row_con_R = []
        for t_dex in range(len(t_con_f)):
            if t_con_f[t_dex] >= next_t:
                to_add_row_con_E.append(E_con_f[t_dex])
                to_add_row_con_I.append(I_con_f[t_dex])
                to_add_row_con_R.append(R_con_f[t_dex])
                next_t += 1
        E_con_full.append(to_add_row_con_E)
        I_con_full.append(to_add_row_con_I)
        R_con_full.append(to_add_row_con_R)
    # Treatment clusters
    E_trt_full = []
    I_trt_full = []
    R_trt_full = []
    for trt_dex in range(n_cluster):
        t_trt_f = sims_trt[chosen_trt[trt_dex]][0]
        E_trt_f = sims_trt[chosen_trt[trt_dex]][2]
        I_trt_f = sims_trt[chosen_trt[trt_dex]][3]
        R_trt_f = sims_trt[chosen_trt[trt_dex]][4]
        next_t = 0
        to_add_row_trt_E = []
        to_add_row_trt_I = []
        to_add_row_trt_R = []
        for t_dex in range(len(t_trt_f)):
            if t_trt_f[t_dex] >= next_t:
                to_add_row_trt_E.append(E_trt_f[t_dex])
                to_add_row_trt_I.append(I_trt_f[t_dex])
                to_add_row_trt_R.append(R_trt_f[t_dex])
                next_t += 1
        E_trt_full.append(to_add_row_trt_E)
        I_trt_full.append(to_add_row_trt_I)
        R_trt_full.append(to_add_row_trt_R)
    # Go through each of the time steps
    EIR_con_sums = []
    EIR_trt_sums = []
    for t in range(200):
        E_con_sum = 0
        E_trt_sum = 0
        I_con_sum = 0
        I_trt_sum = 0
        R_con_sum = 0
        R_trt_sum = 0
        for comm_num in range(n_cluster):
            # Control clusters
            # Get current E
            curr_E_con = 0
            if len(E_con_full[comm_num]) > t:
                curr_E_con = E_con_full[comm_num][t]
            else:
                curr_E_con += E_con_full[comm_num][-1]
            E_con_sum += curr_E_con
            # Get current I
            curr_I_con = 0
            if len(I_con_full[comm_num]) > t:
                curr_I_con = I_con_full[comm_num][t]
            else:
                curr_I_con += I_con_full[comm_num][-1]
            I_con_sum += curr_I_con
            # Get current R
            curr_R_con = 0
            if len(R_con_full[comm_num]) > t:
                curr_R_con = R_con_full[comm_num][t]
            else:
                curr_R_con += R_con_full[comm_num][-1]
            R_con_sum += curr_R_con
            
            # Treatment clusters
            # Get current E
            curr_E_trt = 0
            if len(E_trt_full[comm_num]) > t:
                curr_E_trt = E_trt_full[comm_num][t]
            else:
                curr_E_trt += E_trt_full[comm_num][-1]
            E_trt_sum += curr_E_trt
            # Get current I
            curr_I_trt = 0
            if len(I_trt_full[comm_num]) > t:
                curr_I_trt = I_trt_full[comm_num][t]
            else:
                curr_I_trt += I_trt_full[comm_num][-1]
            I_trt_sum += curr_I_trt
            # Get current R
            curr_R_trt = 0
            if len(R_trt_full[comm_num]) > t:
                curr_R_trt = R_trt_full[comm_num][t]
            else:
                curr_R_trt += R_trt_full[comm_num][-1]
            R_trt_sum += curr_R_trt
        EIR_con_sums.append((E_con_sum + I_con_sum + R_con_sum))
        EIR_trt_sums.append((E_trt_sum + I_trt_sum + R_trt_sum))
    EIR_cons.append(EIR_con_sums)
    EIR_trts.append(EIR_trt_sums)
fig = plt.figure()
ax1 = fig.add_subplot(131)
for i in range(len(EIR_cons)):
    newList = [(x / (n_cluster * 1000)) for x in EIR_cons[i]]
    newList2 = [(x / (n_cluster * 1000)) for x in EIR_trts[i]]
    ax1.plot(list(range(len(EIR_cons[i]))), newList, color='orange', linewidth=0.5, alpha=0.5)
    ax1.plot(list(range(len(EIR_trts[i]))), newList2, color='grey', linewidth=0.5, alpha=0.5)
    ax1.set_title('A')
    ax1.set_ylabel('cum. incidence')
    ax1.set_xlabel('Days since initial infection')
    ax1.text(0.05, 0.95, 'k = 0.4\nβ reduction = 40%', transform=ax1.transAxes, fontsize=17,
        verticalalignment='top', bbox=props)
    ax1.set_xlim([11, 63])
    ax1.set_ylim([0, 0.1])
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
EIR_cons = [] 
EIR_trts = []
n_cluster = 100
for i in range(200):
    chosen_con = random.sample(list(range(len(sims_con))), n_cluster)
    chosen_trt = random.sample(list(range(len(sims_trt))), n_cluster)
    # Control clusters
    E_con_full = []
    I_con_full = []
    R_con_full = []
    for con_dex in range(n_cluster):
        t_con_f = sims_con[chosen_con[con_dex]][0]
        E_con_f = sims_con[chosen_con[con_dex]][2]
        I_con_f = sims_con[chosen_con[con_dex]][3]
        R_con_f = sims_con[chosen_con[con_dex]][4]
        next_t = 0
        to_add_row_con_E = []
        to_add_row_con_I = []
        to_add_row_con_R = []
        for t_dex in range(len(t_con_f)):
            if t_con_f[t_dex] >= next_t:
                to_add_row_con_E.append(E_con_f[t_dex])
                to_add_row_con_I.append(I_con_f[t_dex])
                to_add_row_con_R.append(R_con_f[t_dex])
                next_t += 1
        E_con_full.append(to_add_row_con_E)
        I_con_full.append(to_add_row_con_I)
        R_con_full.append(to_add_row_con_R)
    # Treatment clusters
    E_trt_full = []
    I_trt_full = []
    R_trt_full = []
    for trt_dex in range(n_cluster):
        t_trt_f = sims_trt[chosen_trt[trt_dex]][0]
        E_trt_f = sims_trt[chosen_trt[trt_dex]][2]
        I_trt_f = sims_trt[chosen_trt[trt_dex]][3]
        R_trt_f = sims_trt[chosen_trt[trt_dex]][4]
        next_t = 0
        to_add_row_trt_E = []
        to_add_row_trt_I = []
        to_add_row_trt_R = []
        for t_dex in range(len(t_trt_f)):
            if t_trt_f[t_dex] >= next_t:
                to_add_row_trt_E.append(E_trt_f[t_dex])
                to_add_row_trt_I.append(I_trt_f[t_dex])
                to_add_row_trt_R.append(R_trt_f[t_dex])
                next_t += 1
        E_trt_full.append(to_add_row_trt_E)
        I_trt_full.append(to_add_row_trt_I)
        R_trt_full.append(to_add_row_trt_R)
    # Go through each of the time steps
    EIR_con_sums = []
    EIR_trt_sums = []
    for t in range(200):
        E_con_sum = 0
        E_trt_sum = 0
        I_con_sum = 0
        I_trt_sum = 0
        R_con_sum = 0
        R_trt_sum = 0
        for comm_num in range(n_cluster):
            # Control clusters
            # Get current E
            curr_E_con = 0
            if len(E_con_full[comm_num]) > t:
                curr_E_con = E_con_full[comm_num][t]
            else:
                curr_E_con += E_con_full[comm_num][-1]
            E_con_sum += curr_E_con
            # Get current I
            curr_I_con = 0
            if len(I_con_full[comm_num]) > t:
                curr_I_con = I_con_full[comm_num][t]
            else:
                curr_I_con += I_con_full[comm_num][-1]
            I_con_sum += curr_I_con
            # Get current R
            curr_R_con = 0
            if len(R_con_full[comm_num]) > t:
                curr_R_con = R_con_full[comm_num][t]
            else:
                curr_R_con += R_con_full[comm_num][-1]
            R_con_sum += curr_R_con
            
            # Treatment clusters
            # Get current E
            curr_E_trt = 0
            if len(E_trt_full[comm_num]) > t:
                curr_E_trt = E_trt_full[comm_num][t]
            else:
                curr_E_trt += E_trt_full[comm_num][-1]
            E_trt_sum += curr_E_trt
            # Get current I
            curr_I_trt = 0
            if len(I_trt_full[comm_num]) > t:
                curr_I_trt = I_trt_full[comm_num][t]
            else:
                curr_I_trt += I_trt_full[comm_num][-1]
            I_trt_sum += curr_I_trt
            # Get current R
            curr_R_trt = 0
            if len(R_trt_full[comm_num]) > t:
                curr_R_trt = R_trt_full[comm_num][t]
            else:
                curr_R_trt += R_trt_full[comm_num][-1]
            R_trt_sum += curr_R_trt
        EIR_con_sums.append((E_con_sum + I_con_sum + R_con_sum))
        EIR_trt_sums.append((E_trt_sum + I_trt_sum + R_trt_sum))
    EIR_cons.append(EIR_con_sums)
    EIR_trts.append(EIR_trt_sums)
ax2 = fig.add_subplot(132)
for i in range(len(EIR_cons)):
    newList = [(x / (n_cluster * 1000)) for x in EIR_cons[i]]
    newList2 = [(x / (n_cluster * 1000)) for x in EIR_trts[i]]
    ax2.plot(list(range(len(EIR_cons[i]))), newList, color='orange', linewidth=0.5, alpha=0.5)
    ax2.plot(list(range(len(EIR_trts[i]))), newList2, color='grey', linewidth=0.5, alpha=0.5)
    ax2.set_title('B')
    ax2.set_xlabel('Days since initial infection')
    ax2.text(0.05, 0.95, 'k = 0.4\nβ reduction = 20%', transform=ax2.transAxes, fontsize=17,
        verticalalignment='top', bbox=props)
    ax2.set_xlim([11, 63])
    ax2.set_ylim([0, 0.1])
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
EIR_cons = [] 
EIR_trts = []
n_cluster = 100
for i in range(200):
    chosen_con = random.sample(list(range(len(sims_con))), n_cluster)
    chosen_trt = random.sample(list(range(len(sims_trt))), n_cluster)
    # Control clusters
    E_con_full = []
    I_con_full = []
    R_con_full = []
    for con_dex in range(n_cluster):
        t_con_f = sims_con[chosen_con[con_dex]][0]
        E_con_f = sims_con[chosen_con[con_dex]][2]
        I_con_f = sims_con[chosen_con[con_dex]][3]
        R_con_f = sims_con[chosen_con[con_dex]][4]
        next_t = 0
        to_add_row_con_E = []
        to_add_row_con_I = []
        to_add_row_con_R = []
        for t_dex in range(len(t_con_f)):
            if t_con_f[t_dex] >= next_t:
                to_add_row_con_E.append(E_con_f[t_dex])
                to_add_row_con_I.append(I_con_f[t_dex])
                to_add_row_con_R.append(R_con_f[t_dex])
                next_t += 1
        E_con_full.append(to_add_row_con_E)
        I_con_full.append(to_add_row_con_I)
        R_con_full.append(to_add_row_con_R)
    # Treatment clusters
    E_trt_full = []
    I_trt_full = []
    R_trt_full = []
    for trt_dex in range(n_cluster):
        t_trt_f = sims_trt[chosen_trt[trt_dex]][0]
        E_trt_f = sims_trt[chosen_trt[trt_dex]][2]
        I_trt_f = sims_trt[chosen_trt[trt_dex]][3]
        R_trt_f = sims_trt[chosen_trt[trt_dex]][4]
        next_t = 0
        to_add_row_trt_E = []
        to_add_row_trt_I = []
        to_add_row_trt_R = []
        for t_dex in range(len(t_trt_f)):
            if t_trt_f[t_dex] >= next_t:
                to_add_row_trt_E.append(E_trt_f[t_dex])
                to_add_row_trt_I.append(I_trt_f[t_dex])
                to_add_row_trt_R.append(R_trt_f[t_dex])
                next_t += 1
        E_trt_full.append(to_add_row_trt_E)
        I_trt_full.append(to_add_row_trt_I)
        R_trt_full.append(to_add_row_trt_R)
    # Go through each of the time steps
    EIR_con_sums = []
    EIR_trt_sums = []
    for t in range(200):
        E_con_sum = 0
        E_trt_sum = 0
        I_con_sum = 0
        I_trt_sum = 0
        R_con_sum = 0
        R_trt_sum = 0
        for comm_num in range(n_cluster):
            # Control clusters
            # Get current E
            curr_E_con = 0
            if len(E_con_full[comm_num]) > t:
                curr_E_con = E_con_full[comm_num][t]
            else:
                curr_E_con += E_con_full[comm_num][-1]
            E_con_sum += curr_E_con
            # Get current I
            curr_I_con = 0
            if len(I_con_full[comm_num]) > t:
                curr_I_con = I_con_full[comm_num][t]
            else:
                curr_I_con += I_con_full[comm_num][-1]
            I_con_sum += curr_I_con
            # Get current R
            curr_R_con = 0
            if len(R_con_full[comm_num]) > t:
                curr_R_con = R_con_full[comm_num][t]
            else:
                curr_R_con += R_con_full[comm_num][-1]
            R_con_sum += curr_R_con
            
            # Treatment clusters
            # Get current E
            curr_E_trt = 0
            if len(E_trt_full[comm_num]) > t:
                curr_E_trt = E_trt_full[comm_num][t]
            else:
                curr_E_trt += E_trt_full[comm_num][-1]
            E_trt_sum += curr_E_trt
            # Get current I
            curr_I_trt = 0
            if len(I_trt_full[comm_num]) > t:
                curr_I_trt = I_trt_full[comm_num][t]
            else:
                curr_I_trt += I_trt_full[comm_num][-1]
            I_trt_sum += curr_I_trt
            # Get current R
            curr_R_trt = 0
            if len(R_trt_full[comm_num]) > t:
                curr_R_trt = R_trt_full[comm_num][t]
            else:
                curr_R_trt += R_trt_full[comm_num][-1]
            R_trt_sum += curr_R_trt
        EIR_con_sums.append((E_con_sum + I_con_sum + R_con_sum))
        EIR_trt_sums.append((E_trt_sum + I_trt_sum + R_trt_sum))
    EIR_cons.append(EIR_con_sums)
    EIR_trts.append(EIR_trt_sums)
ax3 = fig.add_subplot(133)
for i in range(len(EIR_cons)):
    newList = [(x / (n_cluster * 1000)) for x in EIR_cons[i]]
    newList2 = [(x / (n_cluster * 1000)) for x in EIR_trts[i]]
    ax3.plot(list(range(len(EIR_cons[i]))), newList, color='orange', linewidth=0.5, alpha=0.5)
    ax3.plot(list(range(len(EIR_trts[i]))), newList2, color='grey', linewidth=0.5, alpha=0.5)
    ax3.set_title('C')
    ax3.set_xlabel('Days since initial infection')
    ax3.text(0.05, 0.95, 'k = 0.1\nβ reduction = 20%', transform=ax3.transAxes, fontsize=17,
        verticalalignment='top', bbox=props)
    ax3.set_xlim([11, 63])
    ax3.set_ylim([0, 0.1])
    ax3.axvline(x=41, color='black', linewidth=0.8, linestyle='-.')
    ax3.axvline(x=30, color='black', linewidth=0.8, linestyle='-.')
    custom_lines = [Line2D([0], [0], color='orange', lw=4),
                    Line2D([0], [0], color='grey', lw=4)]

# Plot the legend -------------------------------------------------------------
lgd = ax2.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2)

# Save figure -----------------------------------------------------------------
fig.savefig(home + '/NPI/code_output/figs/Sheen_Fig1.tiff', dpi=600, bbox_extra_artists=(lgd,), bbox_inches='tight')
fig.savefig(home + '/NPI/code_output/figs/Sheen_Fig1.eps', dpi=600, bbox_extra_artists=(lgd,), bbox_inches='tight')
fig.savefig(home + '/NPI/code_output/figs/Sheen_Fig1.png', dpi=72, bbox_extra_artists=(lgd,), bbox_inches='tight')

        
