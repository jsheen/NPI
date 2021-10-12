#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 20:40:47 2021

@description: script used to create cluster simulations, with and
              without an enacted NPI intervention on day 28 of simulation.

"""
# Import libraries and set seeds ----------------------------------------------
import numpy as np
import random
import matplotlib.pyplot as plt
import statistics
import networkx as nx
from collections import defaultdict, Counter
import EoN
import math
from pathlib import Path
home = str(Path.home())
random.seed(1)
gen = np.random.Generator(np.random.PCG64(1))
# Create parameter sets to run ------------------------------------------------
rts = [1.5, 2]
overdispersions = [0.1, 0.4, 0.7]
clusters = [1000]
effects = [0.2, 0.4]
eits = [0.005]
param_sets = []
for i in rts:
    for j in overdispersions:
        for k in clusters:
            for l in effects:
                for m in eits:
                    param_sets.append([i, j, k, l , m])
clusters = [100]
eits = [0.02]
param_sets = []
for i in rts:
    for j in overdispersions:
        for k in clusters:
            for l in effects:
                for m in eits:
                    param_sets.append([i, j, k, l , m])
# For each parameter set, create 3,000 simulations ----------------------------
for param_set in param_sets:
    tgt_R0 = param_set[0]
    k_overdispersion = param_set[1]
    N_cluster = param_set[2]
    effect = param_set[3]
    expected_It_N = param_set[4]
    mean_degree = 15
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
        raise NameError("Algorithm assumes that the expected It / N strictly < initial infections per cluster / N.")
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
    for i in range(500):
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
    beta = statistics.median(beta_lst)
    # Create graphs for input to Gillespie algorithm --------------------------
    H = nx.DiGraph()
    H.add_node('S')
    H.add_edge('E', 'I', rate = 1 / ave_inc_period, weight_label='expose2infect_weight')
    H.add_edge('I', 'R', rate = 1 / ave_inf_period)
    return_statuses = ('S', 'E', 'I', 'R')
    J = nx.DiGraph()
    J.add_edge(('I', 'S'), ('I', 'E'), rate = beta, weight_label='transmission_weight')
    J_treat = nx.DiGraph()
    J_treat.add_edge(('I', 'S'), ('I', 'E'), rate = ((1 - effect) * beta), weight_label='transmission_weight')
    # Set day of intervention at day 28 ---------------------------------------
    med_t = 28
    # Set threshold value of number of infections at time t -------------------
    threshold = 1
    if N_cluster == 1000:
        threshold = 1
    elif N_cluster == 100:
        threshold = 1
    elif N_cluster == 10000:
        threshold = 1
    # Simulate epidemics with/without treatment of effect reduction in beta ---
    nsim = 3000
    It_It1con_It1trt = []
    sim_ctr = 0
    E_It = []
    while (len(It_It1con_It1trt) < nsim):
        if (len(It_It1con_It1trt) % 100 == 0):
            print(len(It_It1con_It1trt))
        sim_ctr += 1
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
        node_attribute_dict = {node: 1 for node in G.nodes()}
        edge_attribute_dict = {edge: 1 for edge in G.edges()}
        nx.set_node_attributes(G, values=node_attribute_dict, name='expose2infect_weight')
        nx.set_edge_attributes(G, values=edge_attribute_dict, name='transmission_weight')
        IC = defaultdict(lambda: 'S')
        for node in range(initial_infections_per_cluster):
            IC[node] = 'I'
        full_first_half = EoN.Gillespie_simple_contagion(G, H, J, IC, return_statuses, tmax = math.ceil(med_t), return_full_data=True) 
        t_first_half = full_first_half.t()
        S_first_half = full_first_half.S()
        E_first_half = full_first_half.summary()[1]['E']
        I_first_half = full_first_half.I()
        R_first_half = full_first_half.R()
        if I_first_half[-1] >= threshold:
            E_It.append(I_first_half[-1])
            nodes_first_half_final = full_first_half.get_statuses(list(G.nodes()), t_first_half[-1])
            curr_IC = defaultdict(lambda: 'S')
            for node in G.nodes():
                 status = nodes_first_half_final[node]
                 curr_IC[node] = status
            full_second_half_con = EoN.Gillespie_simple_contagion(G, H, J, curr_IC, return_statuses, tmax = np.ceil(one_gen_time) * 3, return_full_data=True)    
            t_second_half_con = full_second_half_con.t()
            S_second_half_con = full_second_half_con.S()
            E_second_half_con = full_second_half_con.summary()[1]['E']
            I_second_half_con = full_second_half_con.I()
            R_second_half_con = full_second_half_con.R()
            full_second_half_trt = EoN.Gillespie_simple_contagion(G, H, J_treat, curr_IC, return_statuses, tmax = np.ceil(one_gen_time) * 3, return_full_data=True)    
            t_second_half_trt = full_second_half_trt.t()
            S_second_half_trt = full_second_half_trt.S()
            E_second_half_trt = full_second_half_trt.summary()[1]['E']
            I_second_half_trt = full_second_half_trt.I()
            R_second_half_trt = full_second_half_trt.R()
            one_gen_S_con, one_gen_E_con, one_gen_I_con, one_gen_R_con = (None,) * 4
            two_gen_S_con, two_gen_E_con, two_gen_I_con, two_gen_R_con = (None,) * 4
            three_gen_S_con, three_gen_E_con, three_gen_I_con, three_gen_R_con = (None,) * 4
            one_gen_S_trt, one_gen_E_trt, one_gen_I_trt, one_gen_R_trt = (None,) * 4
            two_gen_S_trt, two_gen_E_trt, two_gen_I_trt, two_gen_R_trt = (None,) * 4
            three_gen_S_trt, three_gen_E_trt, three_gen_I_trt, three_gen_R_trt = (None,) * 4
            if_con_gen_one = 0
            if_con_gen_two = 0
            for t_dex_con in range(len(t_second_half_con)):
                if t_second_half_con[t_dex_con] >= np.ceil(one_gen_time) and (one_gen_S_con is None):
                    one_gen_S_con = S_second_half_con[t_dex_con]
                    one_gen_E_con = E_second_half_con[t_dex_con]
                    one_gen_I_con = I_second_half_con[t_dex_con]
                    one_gen_R_con = R_second_half_con[t_dex_con]
                    if_con_gen_one += 1
                if t_second_half_con[t_dex_con] >= np.ceil(one_gen_time) * 2 and (two_gen_S_con is None):
                    two_gen_S_con = S_second_half_con[t_dex_con]
                    two_gen_E_con = E_second_half_con[t_dex_con]
                    two_gen_I_con = I_second_half_con[t_dex_con]
                    two_gen_R_con = R_second_half_con[t_dex_con]
                    if_con_gen_two += 1
            if if_con_gen_one > 1 or if_con_gen_two > 1:
                raise NameError("If statement used more than once. (control)")
            three_gen_S_con = S_second_half_con[-1]
            three_gen_E_con = E_second_half_con[-1]
            three_gen_I_con = I_second_half_con[-1]
            three_gen_R_con = R_second_half_con[-1]
            if one_gen_S_con is None: one_gen_S_con = S_second_half_con[-1]
            if one_gen_E_con is None: one_gen_E_con = E_second_half_con[-1]
            if one_gen_I_con is None: one_gen_I_con = I_second_half_con[-1]
            if one_gen_R_con is None: one_gen_R_con = R_second_half_con[-1]
            if two_gen_S_con is None: two_gen_S_con = S_second_half_con[-1]
            if two_gen_E_con is None: two_gen_E_con = E_second_half_con[-1]
            if two_gen_I_con is None: two_gen_I_con = I_second_half_con[-1]
            if two_gen_R_con is None: two_gen_R_con = R_second_half_con[-1]
            if_trt_gen_one = 0
            if_trt_gen_two = 0
            for t_dex_trt in range(len(t_second_half_trt)):
                if t_second_half_trt[t_dex_trt] >= np.ceil(one_gen_time) and (one_gen_S_trt is None):
                    one_gen_S_trt = S_second_half_trt[t_dex_trt]
                    one_gen_E_trt = E_second_half_trt[t_dex_trt]
                    one_gen_I_trt = I_second_half_trt[t_dex_trt]
                    one_gen_R_trt = R_second_half_trt[t_dex_trt]
                    if_trt_gen_one += 1
                if t_second_half_trt[t_dex_trt] >= np.ceil(one_gen_time) * 2 and (two_gen_S_trt is None):
                    two_gen_S_trt = S_second_half_trt[t_dex_trt]
                    two_gen_E_trt = E_second_half_trt[t_dex_trt]
                    two_gen_I_trt = I_second_half_trt[t_dex_trt]
                    two_gen_R_trt = R_second_half_trt[t_dex_trt]
                    if_trt_gen_two += 1
            if if_trt_gen_one > 1 or if_trt_gen_two > 1:
                raise NameError("If statement used more than once. (treatment)")
            three_gen_S_trt = S_second_half_trt[-1]
            three_gen_E_trt = E_second_half_trt[-1]
            three_gen_I_trt = I_second_half_trt[-1]
            three_gen_R_trt = R_second_half_trt[-1]
            if one_gen_S_trt is None: one_gen_S_trt = S_second_half_trt[-1]
            if one_gen_E_trt is None: one_gen_E_trt = E_second_half_trt[-1]
            if one_gen_I_trt is None: one_gen_I_trt = I_second_half_trt[-1]
            if one_gen_R_trt is None: one_gen_R_trt = R_second_half_trt[-1]
            if two_gen_S_trt is None: two_gen_S_trt = S_second_half_trt[-1]
            if two_gen_E_trt is None: two_gen_E_trt = E_second_half_trt[-1]
            if two_gen_I_trt is None: two_gen_I_trt = I_second_half_trt[-1]
            if two_gen_R_trt is None: two_gen_R_trt = R_second_half_trt[-1]
            to_add_row = [S_first_half[-1], E_first_half[-1], I_first_half[-1], R_first_half[-1],
                          one_gen_S_con, one_gen_E_con, one_gen_I_con, one_gen_R_con,
                          two_gen_S_con, two_gen_E_con, two_gen_I_con, two_gen_R_con,
                          three_gen_S_con, three_gen_E_con, three_gen_I_con, three_gen_R_con,
                          one_gen_S_trt, one_gen_E_trt, one_gen_I_trt, one_gen_R_trt,
                          two_gen_S_trt, two_gen_E_trt, two_gen_I_trt, two_gen_R_trt,
                          three_gen_S_trt, three_gen_E_trt, three_gen_I_trt, three_gen_R_trt]
            It_It1con_It1trt.append(to_add_row)
    filename = home + "/NPI/code_output/res/" + str(tgt_R0) + "_" + str(N_cluster) + "_" + str(k_overdispersion) + "_" + str(effect) + "_" + str(expected_It_N) + "_day28.csv"
    with open(filename, 'w') as out_f:
        for sim_dex in range(len(It_It1con_It1trt)):
            for entry_dex in range(len(It_It1con_It1trt[0])):
                out_f.write(str(It_It1con_It1trt[sim_dex][entry_dex]))
                out_f.write(",")
            out_f.write("\n")
        out_f.write(str(sim_ctr))
        out_f.write(",")
        out_f.write(str(statistics.mean(E_It)))
        out_f.write(",")
        out_f.write(str(beta))
        out_f.write(",")
        out_f.write(str(med_t))
