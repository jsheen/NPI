#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  4 17:58:11 2021

@author: Justin
"""
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from collections import defaultdict, Counter
import collections
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
k_overdispersion = 0.1
p = 1.0 - mean_degree / (mean_degree + k_overdispersion)
R0s = []
nccs = []
for i in range(200):
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
                z[i] = 1
        if (sum(z) % 2 == 0):
            continue_loop = False
    G=nx.configuration_model(z)
    G=nx.Graph(G)
    G.remove_edges_from(nx.selfloop_edges(G))
    beta=0.04
    est_R0 = estimate_R0(G, tau=beta, gamma=1/5)
    R0s.append(est_R0)
    ncc = nx.number_connected_components(G)
    nccs.append(ncc / N_cluster)
plt.hist(nccs)
#plt.hist(R0s)
'''
plt.title(str(k_overdispersion))
ccs = nx.connected_components(G)
size_ccs = []
for cc in ccs:
    size_ccs.append(G.subgraph(cc).copy().number_of_nodes())
plt.hist(size_ccs)
degree_sequence = sorted([d for n, d in G.degree()], reverse=True)  # degree sequence
degreeCount = collections.Counter(degree_sequence)
deg, cnt = zip(*degreeCount.items())

fig, ax = plt.subplots()
plt.bar(deg, cnt, width=0.80, color="b")
'''
'''
k_overdispersion = 0.7
p = 1.0 - mean_degree / (mean_degree + k_overdispersion)
test = []
for i in range(1000):
    test.append(np.random.negative_binomial(k_overdispersion, p))
plt.hist(test)
plt.title(str(k_overdispersion))
'''