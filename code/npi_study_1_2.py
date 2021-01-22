# -*- coding: utf-8 -*-
"""NPI_study_1.2.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/18KBbJxn5bjdltnh6Z4mQt5XmPCFqvkRg

# Non-Pharmaceutical Intervention (NPI) cluster randomized study and stepped-wedge analysis randomized study using EoN
- Adaptation of Outbreak_simulations.R code by Lee et al. Code here: https://github.com/leekshaffer/SW-CRT-outbreak
- Primary use is for a speed up in simulation time using Joel C. Miller's Gillespie algorithm implementation in the EoN package. Documentation here: https://epidemicsonnetworks.readthedocs.io/en/latest/index.html
"""

# Load libraries ---------------------------------------------------------------
#!pip install EoN
import EoN
import networkx as nx
from collections import defaultdict
import random
import math
import numpy as np
import os
from os import path
from pathlib import Path

effects = [[1, 100, 0.05, 0.5, 500, 28],
           [1, 100, 0.05, 0.6, 500, 28]]

for effect in effects:
    # Set seed ---------------------------------------------------------------------
    random.seed(0)
    np.random.seed(0)
    
    # Number of simulated trials ---------------------------------------------------
    nsim = 300
    
    # Directory to save results ----------------------------------------------------
    directory_plots = "/Users/jsheen/SW-CRT-outbreak/NPI_study/EoN/code_output/csvs_1_2/"
    
    """# Changeable simulation parameters
    - Population structure parameters, Epidemic parameters, Trial start and end dates, General CRT parameters, Parallel-Arm CRT parameters, SWT parameters
    """
    
    # Population structure parameters ----------------------------------------------
    num_communities = effect[1] # Number of communities
    num_clusters_enrolled_per_day = effect[1] # Num. clusters targeted for enrollment. Must be <= to the number of communities
    ave_community_size = effect[4] # Average size of one community
    community_size_range = 0 # Range of community sizes (sizes are uniformly distributed on this range)
    rate_within = 0.03 # Probability of an edge between two nodes in the same community
    rate_between = 0 # Probability of an edge between two nodes in different communities
    
    # SEIR epidemic parameters -----------------------------------------------------
    direct_NPIE = effect[3] # Leaky multiplicative efficacy of NPI
    beta = effect[2] # Per-time-step hazard of infection for a susceptible nodes from an infectious neighbour
    NPIE_beta = (1 - direct_NPIE) * beta # Per-time-step hazard of infection for a susceptible NPIE node from an infectious neighbour
    incperiod_shape = 5.807 # Gamma-distribution parameters of incubation and 
                            # infectious period. Parameters are from 
                            # Lauer et al. 2020: https://www.acpjournals.org/doi/10.7326/M20-0504
    incperiod_rate = 1 / 0.948
    infperiod_shape = 1.13
    infperiod_rate = 0.226
    ave_inc_period = incperiod_shape / incperiod_rate
    ave_inf_period = infperiod_shape / infperiod_rate
    FixInc = 0
    
    # Starting distribution parameters ---------------------------------------------
    lambda_poisson = 2
    
    # Constant external forcing ----------------------------------------------------
    constant_rate_importation_u = 0
    constant_rate_importation_t = 0
    
    # Calculate R0 -----------------------------------------------------------------
    R0 = (1 - (infperiod_rate / (infperiod_rate+beta)) ** infperiod_shape) * (((ave_community_size - 1) * (1 - rate_within) * rate_within + 
         (num_communities - 1) * ave_community_size * (1 - rate_between) * rate_between + 
         ((ave_community_size-1) * rate_within + (num_communities - 1) * ave_community_size * rate_between) ** 2)/
         ((ave_community_size-1) * rate_within + (num_communities - 1) * ave_community_size * rate_between) - 1)
    
    # Trial start and end date -----------------------------------------------------
    trial_start_day = effect[5] # First day of trial enrollment, relative to start of epidemic

    # General CRT Parameters -------------------------------------------------------
    treatment_cluster_coverage = effect[0]
    
    # Parallel-Arm CRT Parameters --------------------------------------------------
    if (num_clusters_enrolled_per_day > num_communities):
       raise NameError("Enrolling too many communities!")
    enrollment_period = 1 # Num. of days over which subjects are enrolled
    
    # CRT or SWT trial design ------------------------------------------------------
    trial_design = "CRT"
    
    """# Methods to set up randomized-control trial
    
    ## Set SEIR epidemic parameters
    """
    
    # Describe transition rates through H and J graphs for SEIR epidemic -----------
    H = nx.DiGraph()
    H.add_node('S_u')
    H.add_node('S_t')
    H.add_edge('S_u', 'E', rate = constant_rate_importation_u, weight_label='susc_u2expose_weight')
    H.add_edge('S_t', 'E', rate = constant_rate_importation_t, weight_label='susc_t2expose_weight')
    H.add_edge('E', 'I', rate = 1 / ave_inc_period, weight_label='expose2infect_weight')
    H.add_edge('I', 'R', rate = 1 / ave_inf_period, weight_label='infect2recover_weight')
    
    J = nx.DiGraph()
    J.add_edge(('I', 'S_u'), ('I', 'E'), rate = beta, weight_label='transmission_weight_untreated')
    J.add_edge(('I', 'S_t'), ('I', 'E'), rate = NPIE_beta, weight_label='transmission_weight_treated')
    
    return_statuses = ('S_u', 'S_t', 'E', 'I', 'R')
    
    """## Create trial_network
    - Create isolated_network structure
    - Add edges between isolated_network's connected components
    - Create trial_network, which is the isolated_network with edges added between connected components
    """
    
    def create_trial_network():
        # Create isolated communities ----------------------------------------------
        isolated_network = nx.Graph()
        communities = []
        for comm_num in range(num_communities):
            # Randomly select the community size -----------------------------------
            comm_size = np.random.randint(low=ave_community_size - (community_size_range / 2),
                                          high=ave_community_size + (community_size_range / 2) + 1,
                                          size=1)[0]
        
            # Create the community -------------------------------------------------
            comm_network = nx.fast_gnp_random_graph(comm_size, rate_within)
    
            # Check sequential-ness of comm_network --------------------------------
            should_be = 0
            for node in np.sort(list(comm_network.nodes())):
                if node != should_be:
                    print(np.sort(list(comm_network.nodes())))
                    print("Sequential-ness broken!")
                should_be += 1
    
            # Relabeling of nodes required for all communities besides the first ---
            if (len(communities) > 0):
                new_node_names = dict()
                counter = 0
                for old_node_name in list(comm_network.nodes()):
                    new_node_names[old_node_name] = (counter + max(list(isolated_network.nodes())) + 1)
                    counter += 1
                comm_network = nx.relabel_nodes(comm_network, new_node_names)
    
            # Add the community to the isolated_network ----------------------------
            isolated_network.add_nodes_from(comm_network.nodes())
            isolated_network.add_edges_from(comm_network.edges())
    
            # Add the community to the communities list ----------------------------
            communities.append(comm_network)
    
        # Check for duplicates in communities --------------------------------------
        for comm_dex_one in range(len(communities)):
            for comm_dex_two in range((comm_dex_one + 1), len(communities)):
                if len(set(list(communities[comm_dex_one].nodes())).intersection(set(list(communities[comm_dex_two].nodes())))) > 0:
                    print(communities[comm_dex_one].nodes())
                    print(communities[comm_dex_two].nodes())
                    raise NameError("Duplicate in communities: " + str(comm_dex_one) + " and " + str(comm_dex_two))
    
        # Add between community edges if needed ------------------------------------
        trial_network = isolated_network.copy()
        if (rate_between > 0):
            for cc_one in range(len(communities)):
                for cc_two in range(cc_one, len(communities)):
                    for node_cc_one in cc_one.nodes():
                        for node_cc_two in cc_two.nodes():
                            if (np.random.binomial(n=1, p=rate_between, size=1)[0] == 1):
                                trial_network.add_edge_from(node_cc_one, node_cc_two)
    
        # Check for sequential-ness ------------------------------------------------
        should_be = 0
        for node in np.sort(list(trial_network.nodes())):
            if node != should_be:
                print(trial_network.nodes())
                raise NameError("Sequential-ness broken.")
            should_be += 1
    
        return communities, trial_network
    
    """## Set transitions on trial_network
    - Set the spontaneous transition rate parameters and transmission rate parameters for the SEIR simulation
    """
    
    # Method to add transitions ----------------------------------------------------
    def set_transitions(trial_network):
        # Add infectious and incubation periods for each node by randomly drawing
        # from the gamma distribution described by the parameters incperiod_shape, 
        # incperiod_rate, infperiod_shape, and, infperiod_rate ---------------------
        node_attribute_inc_dict = {node: 1 for node in trial_network.nodes()}
        node_attribute_inf_dict = {node: 1 for node in trial_network.nodes()}
        
        # Add spontaneous introductions for each node ------------------------------
        node_attribute_intro_untreat_dict = {node: 1 for node in trial_network.nodes()}
        node_attribute_intro_treat_dict = {node: 1 for node in trial_network.nodes()}
        
        # Set node attributes ------------------------------------------------------
        nx.set_node_attributes(trial_network, values=node_attribute_intro_untreat_dict, name='susc_u2expose_weight')
        nx.set_node_attributes(trial_network, values=node_attribute_intro_treat_dict, name='susc_t2expose_weight')
        nx.set_node_attributes(trial_network, values=node_attribute_inc_dict, name='expose2infect_weight')
        nx.set_node_attributes(trial_network, values=node_attribute_inf_dict, name='infect2recover_weight')
    
        # Add transmission weight to edges for untreated individuals ---------------
        edge_untreat_attribute_dict = {edge: 1 for edge in trial_network.edges()}
        edge_treat_attribute_dict = {edge: 1 for edge in trial_network.edges()}
    
        # Set edge attributes ------------------------------------------------------
        nx.set_edge_attributes(trial_network, values=edge_untreat_attribute_dict, name='transmission_weight_untreated')
        nx.set_edge_attributes(trial_network, values=edge_treat_attribute_dict, name='transmission_weight_treated')
        
        return node_attribute_inc_dict, node_attribute_inf_dict, node_attribute_intro_untreat_dict, node_attribute_intro_treat_dict, edge_untreat_attribute_dict, edge_treat_attribute_dict
    
    """## Pre-assign nodes to treatment or control for a cluster-randomized trial
    - Enroll participants from each community in trial
    - Assign to treatment or control groups
    - Record details of trial in a global dataframe
    """
    
    # Pre-assign communities to either treatment or control ------------------------
    def assign_trial(communities):
        # Assign each community to either treatment or control (i.e. 1 or 0) -------
        treated = random.sample(list(range(len(communities))), math.floor(len(communities) / 2))
        trial_assignment = [-1] * len(communities)
        for community_dex in range(len(trial_assignment)):
            if community_dex in treated:
                trial_assignment[community_dex] = 1
            else:
                trial_assignment[community_dex] = 0
    
        # Check that the number of communities is correct --------------------------
        if len(communities) != num_communities: raise NameError("Wrong num. communities!")
        if len(trial_assignment) != num_communities: raise NameError("Wrong num. communities!")
    
        return trial_assignment
    
    """## Schedule when treatment groups will receive treatment
    - Schedule is dependent on trial design (i.e. CRT or SWT)
    """
    
    def schedule_treatment(trial_assignment):
        treatment_schedule = []
        if (trial_design == "CRT"):
            treatment_schedule.append(trial_start_day)
        elif (trial_design == "SWT"):
            raise NameError("SWT not ready yet to be implemented.")
        return treatment_schedule
    
    """## Enroll susceptible individuals for cluster-randomized trial
    - Enroll "number_enrolled_per_community" susceptible individuals from each cluster
    """
    
    def enroll_participants(communities, susceptibles):
        enrolled_participants = []
        for comm_dex in range(len(communities)):
            susc_in_comm = list(set(communities[comm_dex].nodes()).intersection(set(susceptibles)))
            if len(susc_in_comm) < 2:
                raise NameError("Not enough susceptibles available in community for proper trial.")
            
            else:
                enrolled = susc_in_comm

            enrolled_participants.append(enrolled)
        
        if len(enrolled_participants) != num_communities: raise NameError("Wrong num. communities!")
    
        return enrolled_participants
    
    """# Run trial simulations
    
    ## Loop to run trial simulations
    - Assume constant forcing of infections from external population
    """
    
    if not (path.exists(directory_plots)):
        os.mkdir(directory_plots)
    
    if not (path.exists(directory_plots + "csvs/")):
        os.mkdir(directory_plots + "csvs/")
    
    file_name_params = "beta_" + str(beta) + "/effect_" + str(direct_NPIE) + "/intervention_" + str(trial_start_day)
    batch_folder_name = directory_plots + "csvs/" + file_name_params + "/"
    if not (path.exists(batch_folder_name)):
        Path(batch_folder_name).mkdir(parents=True, exist_ok=True)
    
    cumul_to_save_csv = []
    for sim_num in range(nsim):
        print('Iteration: ' + str(sim_num))
        # Create trial_network -----------------------------------------------------
        communities, trial_network = create_trial_network()
        N = trial_network.number_of_nodes()
    
        # Set transitions of trial_network -----------------------------------------
        node_attribute_inc_dict, node_attribute_inf_dict, node_intro_untreat_dict, node_intro_treat_dict, edge_untreat_attribute_dict, edge_treat_attribute_dict = set_transitions(trial_network)
        
        # Assign clusters to treatment or control ----------------------------------
        trial_assignment = assign_trial(communities)
    
        # Check if there are correct number of assignments -------------------------
        num_to_be_treated = trial_assignment.count(1)
        if (num_to_be_treated != math.floor(len(communities) / 2)):
            raise NameError("Number of clusters assigned to treatment not correct.")
        
        # Schedule when treatment groups receive treatment -------------------------
        treatment_schedule = schedule_treatment(trial_assignment)
        
        # Initialize that everyone is initially susceptible in the study population
        # as well as compartment trajectories --------------------------------------
        curr_IC = {node: 'S_u' for node in trial_network.nodes()}
        t = None
        S_u = None
        S_t = None
        E = None
        I = None
        R = None
    
        # Set initial infections for each community --------------------------------
        for community in communities:
            num_init_infecteds = 0
            while (num_init_infecteds <= 1):
                num_init_infecteds = np.random.poisson(lambda_poisson)
            init_infecteds = random.sample(list(community.nodes()), num_init_infecteds)
            for init_infected in init_infecteds:
                curr_IC[init_infected] = "I"
    
        # Split epidemic in to two: up to treatment day and after treatment day ----
        full_first_half = EoN.Gillespie_simple_contagion(trial_network, H, J, curr_IC, return_statuses, tmax = trial_start_day, return_full_data=True)    
        t_first_half = full_first_half.t()
        S_u_first_half = full_first_half.summary()[1]['S_u']
        S_t_first_half = full_first_half.summary()[1]['S_t']
        E_first_half = full_first_half.summary()[1]['E']
        I_first_half = full_first_half.I()
        R_first_half = full_first_half.R()
        
        # Get currently susceptible individuals ------------------------------------
        nodes_first_half_final = full_first_half.get_statuses(list(trial_network.nodes()), t_first_half[-1])
        curr_IC = defaultdict(lambda: 'S_u')
        for node in trial_network.nodes():
             status = nodes_first_half_final[node]
             curr_IC[node] = status
        susceptibles = []
        for key in curr_IC.keys():
            if curr_IC[key] == 'S_u' or curr_IC[key] == 'S_t':
                susceptibles.append(key)
    
        # Enroll susceptible individuals from each cluster -------------------------
        enrolled_participants = enroll_participants(communities, susceptibles)
    
        # Find communities that need treatment applied -----------------------------
        comm_treat_dexs = [i for i, x in enumerate(trial_assignment) if x == 1]
    
        # Change status from S_u to S_t for each individual of a treatment community
        for comm_treat_dex in comm_treat_dexs:
            to_treat_nodes = random.sample(list(communities[comm_treat_dex].nodes()), math.floor(treatment_cluster_coverage * len(communities[comm_treat_dex].nodes())))
            for to_treat_node in to_treat_nodes:
                # Add assigned treated nodes that are S_u to S_t -------------------
                if curr_IC[to_treat_node] == 'S_u':
                    curr_IC[to_treat_node] = 'S_t'
    
        # Run second half of simulation. Only difference with first call is the
        # modified curr_IC argument ------------------------------------------------
        full_second_half = EoN.Gillespie_simple_contagion(trial_network, H, J, curr_IC, return_statuses, tmax = float('Inf'), return_full_data=True)    
        t_second_half = full_second_half.t()
        S_u_second_half = full_second_half.summary()[1]['S_u']
        S_t_second_half = full_second_half.summary()[1]['S_t']
        E_second_half = full_second_half.summary()[1]['E']
        I_second_half = full_second_half.I()
        R_second_half = full_second_half.R()
    
        # Concatenate first half with second half ----------------------------------
        t = np.concatenate((t_first_half, (t_second_half + t_first_half[-1])), axis=None)
        S_u = np.concatenate((S_u_first_half, S_u_second_half), axis=None)
        S_t = np.concatenate((S_t_first_half, S_t_second_half), axis=None)
        E = np.concatenate((E_first_half, E_second_half), axis=None)
        I = np.concatenate((I_first_half, I_second_half), axis=None)
        R = np.concatenate((R_first_half, R_second_half), axis=None)
        
        # Initialize dictionaries of control and treatment time series by community -
        control_t_series = dict()
        treatment_t_series = dict()
        for comm_dex in range(len(communities)):
            treatment_or_not = trial_assignment[comm_dex]
            if treatment_or_not == 1:
                treatment_t_series[comm_dex] = list()
            else:
                control_t_series[comm_dex] = list()
        
        # Loop through first half of simulations -------------------------------------
        ts_first_half = list(range(int(np.floor(max(t_first_half) + 1))))
        ts_second_half = list(range(int(np.floor(max(t_second_half) + 1))))
        for curr_t in ts_first_half:
            status_dict = full_first_half.get_statuses(list(trial_network.nodes()), curr_t)
            for comm_dex in range(len(communities)):
                treatment_or_not = trial_assignment[comm_dex]
                curr_t_S = 0
                curr_t_E = 0
                curr_t_I = 0
                curr_t_R = 0
                for node in communities[comm_dex].nodes():
                    status_of_node = status_dict[node]
                    if status_of_node == "I":
                        curr_t_I += 1
                    elif status_of_node == "E":
                        curr_t_E += 1
                    elif status_of_node == "R":
                        curr_t_R += 1
                    else:
                        curr_t_S += 1
                to_add_cell = str(curr_t) + "_" + str(curr_t_S) + "_" + str(curr_t_E) + "_" + str(curr_t_I) + "_" + str(curr_t_R)
                if treatment_or_not == 1:
                    treatment_t_series[comm_dex].append(to_add_cell)
                else:
                    control_t_series[comm_dex].append(to_add_cell)
        
        # Loop through second half of simulations -------------------------------------
        for curr_t in ts_second_half:
            status_dict = full_second_half.get_statuses(list(trial_network.nodes()), curr_t)
            for comm_dex in range(len(communities)):
                treatment_or_not = trial_assignment[comm_dex]
                curr_t_S = 0
                curr_t_E = 0
                curr_t_I = 0
                curr_t_R = 0
                for node in communities[comm_dex].nodes():
                    status_of_node = status_dict[node]
                    if status_of_node == "I":
                        curr_t_I += 1
                    elif status_of_node == "E":
                        curr_t_E += 1
                    elif status_of_node == "R":
                        curr_t_R += 1
                    else:
                        curr_t_S += 1
                to_add_cell = str(curr_t + ts_first_half[-1] + 1) + "_" + str(curr_t_S) + "_" + str(curr_t_E) + "_" + str(curr_t_I) + "_" + str(curr_t_R)
                if treatment_or_not == 1:
                    treatment_t_series[comm_dex].append(to_add_cell)
                else:
                    control_t_series[comm_dex].append(to_add_cell)
        
        # Create .csv to save of simulation ----------------------------------------
        cumul_to_save_csv.append([control_t_series, treatment_t_series])
    
        # Check that the number of nodes in network is within our expectation ------
        if (trial_network.number_of_nodes() < ((ave_community_size - (community_size_range / 2)) * num_communities) or
            trial_network.number_of_nodes() > ((ave_community_size + (community_size_range / 2) + 1) * num_communities)):
            raise NameError("Number of nodes in network is not within expectation.")
    
        # Check that trial_network has not changed within loop ---------------------
        if N != trial_network.number_of_nodes():
            raise NameError("Number of nodes in trial_network not correct.")
    
        # Check that if rate_between == 0, there are no edges between ccs ----------
        if (rate_between == 0):
            tot_edges_communities = 0
            for community in communities:
                tot_edges_communities += community.number_of_edges()
            if tot_edges_communities != trial_network.number_of_edges():
                raise NameError("Number of edges between communities and trial network are not same.")
    
    # Write nsim number of folders, each having ncomm number of .csvs --------------
    for sim_num in range(len(cumul_to_save_csv)):
        control_dict = cumul_to_save_csv[sim_num][0]
        treatment_dict = cumul_to_save_csv[sim_num][1]
        for comm_num in list(control_dict.keys()):
            csv_name = batch_folder_name + "sim_" + str(sim_num) + "/" 
            Path(csv_name).mkdir(parents=True, exist_ok=True)
            with open(csv_name + "comm_" + str(comm_num) + "_control.csv", 'w') as out_f:
                comm_c = control_dict[comm_num]
                for l in range(len(comm_c)):
                    split_string = comm_c[l].split("_")
                    for l2 in range(len(split_string)):
                        out_f.write(split_string[l2] + ",")
                    out_f.write("\n")
        for comm_num in list(treatment_dict.keys()):
            csv_name = batch_folder_name + "sim_" + str(sim_num) + "/" 
            Path(csv_name).mkdir(parents=True, exist_ok=True)
            with open(csv_name + "comm_" + str(comm_num) + "_treatment.csv", 'w') as out_f:
                comm_t = treatment_dict[comm_num]
                for l in range(len(comm_t)):
                    split_string = comm_t[l].split("_")
                    for l2 in range(len(split_string)):
                        out_f.write(split_string[l2] + ",")
                    out_f.write("\n")
    
    """# Additional tests
    - Testing that networkx can handle two node attributes on the same node
    - Check that rates of transmission and transitions make sense
    - Check that rates of transmissions and transitions reflect heterogeneity
    """
    
    # Testing that networkx can handle two node attributes on the same node --------
    test = nx.fast_gnp_random_graph(2, 1)
    node_one_dict = {node: "test1" for node in test.nodes()}
    node_two_dict = {node: "test2" for node in test.nodes()}
    nx.set_node_attributes(test, values=node_one_dict, name='expose2infect_weight')
    nx.set_node_attributes(test, values=node_two_dict, name='infect2recover_weight')
    #print(nx.get_node_attributes(test, name="expose2infect_weight"))
    #print(nx.get_node_attributes(test, name="infect2recover_weight"))
    
    # Check that rates of transmissions and transitions make sense qualitatively ---
    #print(nx.get_node_attributes(trial_network, name="susc_u2expose_weight")) # No heterogeneity
    #print(nx.get_node_attributes(trial_network, name="susc_t2expose_weight")) # No heterogeneity
    #print(nx.get_node_attributes(trial_network, name="expose2infect_weight")) # Yes heterogeneity
    #print(nx.get_node_attributes(trial_network, name="infect2recover_weight")) # Yes heterogeneity
    #print(nx.get_edge_attributes(trial_network, name="transmission_weight_untreated")) # No heterogeneity
    #print(nx.get_edge_attributes(trial_network, name="transmission_weight_treated")) # No heterogeneity
    
    # Check that rates of transmissions and transitions reflect heterogeneity ------
    if len(np.unique(list(nx.get_node_attributes(trial_network, name="susc_u2expose_weight").values()))) > 1:
        raise NameError("Heterogeneity where there should be no heterogeneity.")
    if len(np.unique(list(nx.get_node_attributes(trial_network, name="susc_t2expose_weight").values()))) > 1:
        raise NameError("Heterogeneity where there should be no heterogeneity.")
    if len(np.unique(list(nx.get_edge_attributes(trial_network, name="transmission_weight_untreated").values()))) > 1:
        raise NameError("Heterogeneity where there should be no heterogeneity.")
    if len(np.unique(list(nx.get_edge_attributes(trial_network, name="transmission_weight_treated").values()))) > 1:
        raise NameError("Heterogeneity where there should be no heterogeneity.")
    
    from itertools import chain

    # Check for connected components matching with communities ---------------------
    ccs = list((trial_network.subgraph(c) for c in nx.connected_components(trial_network)))
    community_nodes = []
    for community in communities:
        community_nodes.append(list(community.nodes()))
    community_nodes = list(chain.from_iterable(community_nodes))
    
    ccs_nodes = []
    for cc in ccs:
        ccs_nodes.append(list(cc.nodes()))
    ccs_nodes = list(chain.from_iterable(ccs_nodes))
    
    counter = 0
    for community in communities:
        if set(list(ccs[counter].nodes())) != set(list(community.nodes())):
            print(counter)
        counter += 1