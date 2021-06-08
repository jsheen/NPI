#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  4 23:27:31 2021

@author: Justin
"""
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import statistics
from collections import defaultdict, Counter
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
N_cluster = 1000
mean_degree = 15
k_overdispersion = 0.7
tgt_R0 = 2
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
        for i in range(len(z)):
            if (z[i] == 0):
                z[i] == 1
        if (sum(z) % 2 == 0):
            continue_loop = False
    G=nx.configuration_model(z)
    G=nx.Graph(G)
    G.remove_edges_from(nx.selfloop_edges(G))
    est_R0=3.3
    beta=0.04
    while est_R0 > tgt_R0:
        beta = beta - 0.0001
        est_R0 = estimate_R0(G, tau=beta, gamma=1/5)
    print(beta)
    print(est_R0)
    beta_lst.append(beta)
plt.hist(beta_lst)
plt.title(str(k_overdispersion))