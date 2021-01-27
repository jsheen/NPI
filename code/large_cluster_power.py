#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 20:40:47 2021

@author: jsheen
"""

import numpy as np
import matplotlib.pyplot as plt
import statistics
import networkx as nx
from collections import defaultdict, Counter
import EoN
import random
import math
random.seed(0)
gen = np.random.Generator(np.random.PCG64(0))

# Global parameters -----------------------------------------------------------
effect = 0.4
expected_It_N = 0.01
N_cluster = 10000
mean_degree = 15
k_overdispersion = 0.5
tgt_R0 = 2.00
initial_infections_per_cluster = 30
incperiod_shape = 5.807
incperiod_rate = 1 / 0.948
infperiod_shape = 1.13
infperiod_rate = 0.226
ave_inc_period = incperiod_shape / incperiod_rate
ave_inf_period = infperiod_shape / infperiod_rate
if expected_It_N <= initial_infections_per_cluster / N_cluster:
    raise NameError("Algorithm assumes that the expected It / N strictly < initial infections per cluster / N.")

# Joel's methods to estimate_R0 without bug -----------------------------------
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

# Pre-analysis of degree distribution -----------------------------------------
p = 1.0 - mean_degree / (mean_degree + k_overdispersion)
res = []
for i in range(N_cluster):
    f_res = 0
    f_res = np.random.negative_binomial(k_overdispersion, p)
    res.append(f_res)
for i in range(len(res)):
    if (res[i] == 0):
        res[i] == 1
plt.hist(res)

# Find median beta ------------------------------------------------------------
# 2.00: 0.0094
beta_lst = []
for i in range(100):
    if (i % 10 == 0):
        print(i)
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
    
    ave_inf_period = 5
    est_R0 = 3.3
    beta = 0.03
    while est_R0 > tgt_R0:
        beta = beta - 0.0001
        est_R0 = estimate_R0(G, tau=beta, gamma=1/ave_inf_period)
    beta_lst.append(beta)
plt.hist(beta_lst)
print("Median beta value for tgt_R0: " + statistics.median(beta_lst))
beta = statistics.median(beta_lst)

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

# Find day on average when there is 1% of active infections -------------------
nsim = 1000
I_series = []
for sim_i in range(nsim):
    if sim_i % 100 == 0:
        print(sim_i)
    p = 1.0 - mean_degree / (mean_degree + k_overdispersion)
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
    t, S, E, I, R = EoN.Gillespie_simple_contagion(G, H, J, IC, return_statuses, tmax = 100)
    next_t = 0
    to_add_row = []
    for t_dex in range(len(t)):
        if t[t_dex] >= next_t:
            to_add_row.append(I[t_dex])
            next_t += 1
    I_series.append(to_add_row)
med_t_one_pct = None
for day_dex in range(len(I_series[0])):
    focal_dist = []
    for I_series_dex in range(len(I_series)):
        if len(I_series[I_series_dex]) > day_dex:
            focal_dist.append(I_series[I_series_dex][day_dex] / N_cluster)
    if statistics.mean(focal_dist) >= expected_It_N:
        med_t_one_pct = day_dex
        break

# Simulate epidemics with and without treatment of effect reduction in beta ---
nsim = 1000
It_It1con_It1trt = []
for sim_i in range(nsim):
    if sim_i % 100 == 0:
        print(sim_i)
    p = 1.0 - mean_degree / (mean_degree + k_overdispersion)
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
    if I_first_half[-1] >= 10:
        nodes_first_half_final = full_first_half.get_statuses(list(G.nodes()), t_first_half[-1])
        curr_IC = defaultdict(lambda: 'S')
        for node in G.nodes():
             status = nodes_first_half_final[node]
             curr_IC[node] = status
        full_second_half_con = EoN.Gillespie_simple_contagion(G, H, J, curr_IC, return_statuses, tmax = 11.5, return_full_data=True)    
        t_second_half_con = full_second_half_con.t()
        S_second_half_con = full_second_half_con.S()
        E_second_half_con = full_second_half_con.summary()[1]['E']
        I_second_half_con = full_second_half_con.I()
        R_second_half_con = full_second_half_con.R()
        full_second_half_trt = EoN.Gillespie_simple_contagion(G, H, J_treat, curr_IC, return_statuses, tmax = 11.5, return_full_data=True)    
        t_second_half_trt = full_second_half_trt.t()
        S_second_half_trt = full_second_half_trt.S()
        E_second_half_trt = full_second_half_trt.summary()[1]['E']
        I_second_half_trt = full_second_half_trt.I()
        R_second_half_trt = full_second_half_trt.R()
        to_add_row = [S_first_half[-1], E_first_half[-1], I_first_half[-1], R_first_half[-1],
                      S_second_half_con[-1], E_second_half_con[-1], I_second_half_con[-1], R_second_half_con[-1],
                      S_second_half_trt[-1], E_second_half_trt[-1], I_second_half_trt[-1], R_second_half_trt[-1]]
        It_It1con_It1trt.append(to_add_row)

# Sample ncluster for each arm ------------------------------------------------
nsim = 10000
ncluster = 63
log_It1It_trt = []
log_It1It_con = []
for sim_i in range(nsim):
    chosen = random.sample(list(range(len(It_It1con_It1trt))), ncluster * 2)
    sum_It_con = 0
    sum_It1_con = 0
    for con_dex in range(ncluster):
        St = It_It1con_It1trt[chosen[con_dex]][0]
        Et = It_It1con_It1trt[chosen[con_dex]][1]
        It = It_It1con_It1trt[chosen[con_dex]][2]
        Rt = It_It1con_It1trt[chosen[con_dex]][3]
        St1 = It_It1con_It1trt[chosen[con_dex]][4]
        Et1 = It_It1con_It1trt[chosen[con_dex]][5]
        It1 = It_It1con_It1trt[chosen[con_dex]][6]
        Rt1 = It_It1con_It1trt[chosen[con_dex]][7]
        sampled_It_con = gen.multivariate_hypergeometric([St, Et, It, Rt], nsample=100, size=1)[0][2]
        sampled_It1_con = gen.multivariate_hypergeometric([St1, Et1, It1, Rt1], nsample=100, size=1)[0][2]
        sum_It_con += sampled_It_con
        sum_It1_con += sampled_It1_con
    sum_It_trt = 0
    sum_It1_trt = 0
    for trt_dex in range(ncluster, ncluster * 2):
        St = It_It1con_It1trt[chosen[trt_dex]][0]
        Et = It_It1con_It1trt[chosen[trt_dex]][1]
        It = It_It1con_It1trt[chosen[trt_dex]][2]
        Rt = It_It1con_It1trt[chosen[trt_dex]][3]
        St1 = It_It1con_It1trt[chosen[trt_dex]][8]
        Et1 = It_It1con_It1trt[chosen[trt_dex]][9]
        It1 = It_It1con_It1trt[chosen[trt_dex]][10]
        Rt1 = It_It1con_It1trt[chosen[trt_dex]][11]
        sampled_It_trt = gen.multivariate_hypergeometric([St, Et, It, Rt], nsample=100, size=1)[0][2]
        sampled_It1_trt = gen.multivariate_hypergeometric([St1, Et1, It1, Rt1], nsample=100, size=1)[0][2]
        sum_It_trt += sampled_It_trt
        sum_It1_trt += sampled_It1_trt
    log_It1It_trt.append(np.log((sum_It1_trt + 1) / (sum_It_trt + 1)))
    log_It1It_con.append(np.log((sum_It1_con + 1) / (sum_It_con + 1)))
    
# Plot difference -------------------------------------------------------------
log_It1It_trt.sort()
log_It1It_con.sort()
p_val_x = log_It1It_con[round(0.05 * len(log_It1It_con))]
power = len(np.where(log_It1It_trt <= p_val_x)[0]) / len(log_It1It_trt)
plt.hist(log_It1It_trt, bins=20)
plt.hist(log_It1It_con, bins=20)
plt.xlabel("log((I_t+1 + 1) / (I_t + 1))")
power


