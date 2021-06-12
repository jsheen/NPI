#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 31 17:15:24 2021

@author: Justin Sheen

@description: script used to create the figures showing the simulations 
              trajectories.

""" 
# Import libraries, set seeds and parameters -------------------------------
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 18})
import networkx as nx
from collections import defaultdict, Counter
import EoN
import random
import math
from matplotlib.lines import Line2D
import pickle
from pathlib import Path
home = str(Path.home())
random.seed(1)
gen = np.random.Generator(np.random.PCG64(1))
nsim = 500
k_overdispersion = 0.4
N_cluster = 1000
effect = 0.4
mean_degree = 15
p = 1.0 - mean_degree / (mean_degree + k_overdispersion)
threshold = 1
initial_infections_per_cluster = 4
beta = 0.0072
med_t_one_pct = 30
incperiod_shape = 5.807
incperiod_rate = 1 / 0.948
infperiod_shape = 1.13
infperiod_rate = 0.226
ave_inc_period = incperiod_shape / incperiod_rate
ave_inf_period = infperiod_shape / infperiod_rate
one_gen_time = ave_inc_period + ave_inf_period
# Joel C Miller's methods to estimate_R0 ----------------------------------
def get_Pk(G):
    Nk = Counter(dict(G.degree()).values())
    Pk = {x:Nk[x]/float(G.order()) for x in Nk.keys()}
    return Pk
def get_PGFPrime(Pk):
    maxk = max(Pk.keys())
    ks = np.linspace(0,maxk, maxk+1)
    Pkarray = np.array([Pk.get(k,0) for k in ks])
    return lambda x: Pkarray.dot(ks*x**(ks-1))
def get_PGFDPrime(Pk):
    maxk = max(Pk.keys())
    ks = np.linspace(0,maxk, maxk+1)
    Pkarray = np.array([Pk.get(k,0) for k in ks])
    return lambda x: Pkarray.dot(ks*(ks-1)*x**(ks-2))
def estimate_R0(G, tau = None, gamma = None):
    transmissibility = tau/(tau+gamma)
    Pk = get_Pk(G)
    psiDPrime = get_PGFDPrime(Pk)
    psiPrime = get_PGFPrime(Pk)
    return transmissibility * psiDPrime(1.)/psiPrime(1.)
# Create graphs for input to Gillespie algorithm ------------------------------
H = nx.DiGraph()
H.add_node('S')
H.add_edge('E', 'I', rate = 1 / ave_inc_period, weight_label='expose2infect_weight')
H.add_edge('I', 'R', rate = 1 / ave_inf_period)
return_statuses = ('S', 'E', 'I', 'R')
J = nx.DiGraph()
J.add_edge(('I', 'S'), ('I', 'E'), rate = beta, weight_label='transmission_weight')
J_treat = nx.DiGraph()
J_treat.add_edge(('I', 'S'), ('I', 'E'), rate = ((1 - effect) * beta), weight_label='transmission_weight')
# Run simulations -------------------------------------------------------------
nsim = 1000
sims_con = []
sims_trt = []
sim_ctr = 0
while (len(sims_con) < nsim):
    if (len(sims_con) % 100 == 0):
        print(len(sims_con))
    sim_ctr += 1
    continue_loop = True
    while (continue_loop):
        z = []
        for i in range(N_cluster):
            deg = 0
            deg = np.random.negative_binomial(k_overdispersion, p)
            z.append(deg)
        for i in range(len(z)):
            if (z[i] == 0):
                z[i] == 1
        if (sum(z) % 2 == 0):
            continue_loop = False
    G=nx.configuration_model(z)
    G=nx.Graph(G)
    G.remove_edges_from(nx.selfloop_edges(G))
    node_attribute_dict = {node: 1 for node in G.nodes()}
    edge_attribute_dict = {edge: 1 for edge in G.edges()}
    nx.set_node_attributes(G, values=node_attribute_dict, name='expose2infect_weight')
    nx.set_edge_attributes(G, values=edge_attribute_dict, name='transmission_weight')
    IC = defaultdict(lambda: 'S')
    for node in range(initial_infections_per_cluster):
        IC[node] = 'I'
    full_first_half = EoN.Gillespie_simple_contagion(G, H, J, IC, return_statuses, tmax = math.ceil(med_t_one_pct), return_full_data=True) 
    t_first_half = full_first_half.t()
    S_first_half = full_first_half.S()
    E_first_half = full_first_half.summary()[1]['E']
    I_first_half = full_first_half.I()
    R_first_half = full_first_half.R()
    if I_first_half[-1] >= threshold:
        nodes_first_half_final = full_first_half.get_statuses(list(G.nodes()), t_first_half[-1])
        curr_IC = defaultdict(lambda: 'S')
        for node in G.nodes():
             status = nodes_first_half_final[node]
             curr_IC[node] = status
        full_second_half_con = EoN.Gillespie_simple_contagion(G, H, J, curr_IC, return_statuses, tmax = float('Inf'), return_full_data=True)    
        t_second_half_con = full_second_half_con.t()
        S_second_half_con = full_second_half_con.S()
        E_second_half_con = full_second_half_con.summary()[1]['E']
        I_second_half_con = full_second_half_con.I()
        R_second_half_con = full_second_half_con.R()
        full_second_half_trt = EoN.Gillespie_simple_contagion(G, H, J_treat, curr_IC, return_statuses, tmax=float('Inf'), return_full_data=True)    
        t_second_half_trt = full_second_half_trt.t()
        S_second_half_trt = full_second_half_trt.S()
        E_second_half_trt = full_second_half_trt.summary()[1]['E']
        I_second_half_trt = full_second_half_trt.I()
        R_second_half_trt = full_second_half_trt.R()
        # Concatenate first half with second half ----------------------------------
        t_con = np.concatenate((t_first_half, (t_second_half_con + t_first_half[-1])), axis=None)
        S_con = np.concatenate((S_first_half, S_second_half_con), axis=None)
        E_con = np.concatenate((E_first_half, E_second_half_con), axis=None)
        I_con = np.concatenate((I_first_half, I_second_half_con), axis=None)
        R_con = np.concatenate((R_first_half, R_second_half_con), axis=None)
        t_trt = np.concatenate((t_first_half, (t_second_half_trt + t_first_half[-1])), axis=None)
        S_trt = np.concatenate((S_first_half, S_second_half_trt), axis=None)
        E_trt = np.concatenate((E_first_half, E_second_half_trt), axis=None)
        I_trt = np.concatenate((I_first_half, I_second_half_trt), axis=None)
        R_trt = np.concatenate((R_first_half, R_second_half_trt), axis=None)
    else:
        t_con = None
        t_trt = None
        S_con = None
        S_trt = None
        E_con = None
        E_trt = None
        I_con = None
        I_trt = None
        R_con = None
        R_trt = None
    sims_con.append([t_con, S_con, E_con, I_con, R_con])
    sims_trt.append([t_trt, S_trt, E_trt, I_trt, R_trt])
with open(home + '/NPI/code_output/figs/sims_con.pickle', 'wb') as handle:
    pickle.dump(sims_con, handle, protocol=pickle.HIGHEST_PROTOCOL)
with open(home + '/NPI/code_output/figs/sims_trt.pickle', 'wb') as handle:
    pickle.dump(sims_trt, handle, protocol=pickle.HIGHEST_PROTOCOL)

# Plotting --------------------------------------------------------------------
with open(home + '/NPI/code_output/figs/sims_con.pickle', 'rb') as handle:
    sims_con = pickle.load(handle)
with open(home + '/NPI/code_output/figs/sims_trt.pickle', 'rb') as handle:
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
I_cons = [] # Create 100 simulations, randomly picking n_cluster number of clusters per arm
I_trts = []
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
    I_cons.append(I_con_sums)
    I_trts.append(I_trt_sums)
fig = plt.figure()
ax1 = fig.add_subplot(111)
for i in range(len(I_cons)):
    newList = [(x / (n_cluster * 1000)) for x in I_cons[i]]
    newList2 = [(x / (n_cluster * 1000)) for x in I_trts[i]]
    ax1.plot(list(range(len(I_cons[i]))), newList, color='orange', linewidth=0.5, alpha=0.5)
    ax1.plot(list(range(len(I_trts[i]))), newList2, color='grey', linewidth=0.5, alpha=0.5)
    ax1.set_ylabel('Prevalence')
    ax1.set_xlim([10.5, 100])
    ax1.set_ylim([0, 0.01])
    plt.axvline(x=41, color='black', linewidth=0.8, linestyle='-.')
    plt.axvline(x=30, color='black', linewidth=0.8, linestyle='-.')
    custom_lines = [Line2D([0], [0], color='orange', lw=4),
                    Line2D([0], [0], color='grey', lw=4)]
plt.tight_layout()
plt.savefig(home + '/NPI/code_output/figs/sim_fig_top.png', dpi=300)

        
