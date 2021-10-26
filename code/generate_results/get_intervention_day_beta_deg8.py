#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 19:19:27 2021

@author: Justin Sheen

@description: sensitivity analysis
"""
# Import libraries and set seeds ----------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import statistics
import networkx as nx
from collections import defaultdict, Counter
import EoN
from pathlib import Path
home = str(Path.home())
# Create parameter sets to run ------------------------------------------------
rts = [2]
overdispersions = [0.7]
"""
N_clusters = [10000]
eits = [0.005]
param_sets = []
for i in rts:
    for j in overdispersions:
        for k in N_clusters:
            for l in eits:
               	if (j == 0.1 and k == 10000):
                       param_sets.append([i, j, k, 0.0045]) 
                else:
                       param_sets.append([i, j, k, l])
"""
N_clusters = [100]
eits = [0.02]
param_sets = []
for i in rts:
    for j in overdispersions:
        for k in N_clusters:
            for l in eits:
                param_sets.append([i, j, k, l])

# For each parameter set, create 3,000 simulations ----------------------------
for param_set in param_sets:
    tgt_R0 = param_set[0]
    k_overdispersion = param_set[1]
    N_cluster = param_set[2]
    expected_It_N = param_set[3]
    mean_degree = 8
    initial_infections_per_cluster = None
    if N_cluster == 100:
        initial_infections_per_cluster = 1
    elif N_cluster == 1000:
        initial_infections_per_cluster = 4
    else:
        initial_infections_per_cluster = 40
    incperiod_shape = 5.807
    incperiod_rate = 1 / 0.948
    infperiod_shape = 1.13
    infperiod_rate = 0.226
    ave_inc_period = incperiod_shape / incperiod_rate
    ave_inf_period = infperiod_shape / infperiod_rate
    one_gen_time = ave_inc_period + ave_inf_period
    if expected_It_N <= initial_infections_per_cluster / N_cluster:
        raise NameError("Script assumes expected It / N strictly < initial infections per cluster / N.")
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
    # Find median beta that leads to desired tgt_R0 (500 sims) ----------------
    p = 1.0 - mean_degree / (mean_degree + k_overdispersion)
    beta_lst = []
    for i in range(2000):
        if (i % 100 == 0):
            print(i)
        continue_loop = True
        while (continue_loop):
            z = []
            for i in range(N_cluster):
                deg = 0
                deg = np.random.negative_binomial(k_overdispersion, p)
                z.append(deg)
            if (sum(z) % 2 == 0):
                continue_loop = False
        G=nx.configuration_model(z)
        G=nx.Graph(G)
        G.remove_edges_from(nx.selfloop_edges(G))
        est_R0=3.3
        beta=0.04
        while est_R0 > tgt_R0:
            beta = beta - 0.0001
            est_R0 = estimate_R0(G, tau=beta, gamma=1/ave_inf_period)
        beta_lst.append(beta)
    plt.hist(beta_lst)
    print("Median beta value for tgt_R0: " + str(statistics.median(beta_lst)))
    final_beta = statistics.median(beta_lst)
    # Create graphs for input to Gillespie algorithm --------------------------
    H = nx.DiGraph()
    H.add_node('S')
    H.add_edge('E', 'I', rate = 1 / ave_inc_period, weight_label='expose2infect_weight')
    H.add_edge('I', 'R', rate = 1 / ave_inf_period)
    return_statuses = ('S', 'E', 'I', 'R')
    J = nx.DiGraph()
    J.add_edge(('I', 'S'), ('I', 'E'), rate = final_beta, weight_label='transmission_weight')
    # Find day on average when expected_It_N of active infections (2000 sims) -
    nsim = 2000
    I_series = []
    while (len(I_series) < nsim):
        if (len(I_series) % 200 == 0):
            print(len(I_series))
        continue_loop = True
        while (continue_loop):
            z = []
            for i in range(N_cluster):
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
        t, S, E, I, R = EoN.Gillespie_simple_contagion(G, H, J, IC, return_statuses, tmax = 200)
        next_t = 0
        to_add_row = []
        for t_dex in range(len(t)):
            if t[t_dex] >= next_t:
                to_add_row.append(I[t_dex])
                next_t += 1
        I_series.append(to_add_row)
    interrupt_t = None
    # Find first day of sim where the ave. num. of infects >= expected_It_N ---
    for day_dex in range(nsim):
        focal_dist = []
        for I_series_dex in range(len(I_series)):
            if len(I_series[I_series_dex]) > day_dex:
                focal_dist.append(I_series[I_series_dex][day_dex] / N_cluster)
        if len(focal_dist) <= 200:
            raise NameError("Not enough simulations (<10%) to get average number of infections on this day.")
        print(len(focal_dist))
        print(statistics.mean(focal_dist))
        if statistics.mean(focal_dist) >= expected_It_N:
            interrupt_t = day_dex
            break
    # Write output ------------------------------------------------------------    
    filename = home + "/NPI/code_output/prelim/" + str(tgt_R0) + "_" + str(N_cluster) + "_" + str(k_overdispersion) + "_" + str(expected_It_N) + "_deg8.csv"
    with open(filename, 'w') as out_f:
        out_f.write(str(final_beta))
        out_f.write(",")
        out_f.write(str(interrupt_t))
