#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 23:47:59 2021
@description: script used to create the figures showing the simulations 
              trajectories for lower effect size.
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
nsim = 1000
k_overdispersion = 0.4
N_cluster = 1000
effect = 0.2
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
with open(home + '/NPI/code_output/figs/Sheen_Fig1_pickles/sims_con_lowE.pickle', 'wb') as handle:
    pickle.dump(sims_con, handle, protocol=pickle.HIGHEST_PROTOCOL)
with open(home + '/NPI/code_output/figs/Sheen_Fig1_pickles/sims_trt_lowE.pickle', 'wb') as handle:
    pickle.dump(sims_trt, handle, protocol=pickle.HIGHEST_PROTOCOL)
        
