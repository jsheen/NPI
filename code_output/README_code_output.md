# Output of simulation code

Folder 1: res (simulation results of paper)
- Naming convention of simulation results: R0_N_k_effect_E[I] (_day28 if day28 file, _deg8 if deg8 file, _deg30 if deg30 file)
- Naming convention of search algorithm results: res_R0_N_k_effect_E[I]_nsample (one of eight types of files: ttest_rerun, ttest results; ttest_gen2, ttest results when sampling two generations after intervention; ttest_gen3, ttest results when sampling three generations after intervention; ttest_day28, ttest results when intervention occurs on day 28 of epidemic; M_S, ttest results when matched on number of sampled susceptible individuals at time of intervention; M_SER, ttest results when matched on number of sampled non-infectious individuals at time of intervention; ttest_deg8, ttest results when average degree is 8; ttest_deg30, ttest results when average degree is 30)

Folder 2: tables (tables of paper)
- Master tables of results. Eight types: ttest, ttest results; gen2, ttest_gen2 results when sampling two generations after intervention; ttest_gen3, ttest results when sampling three generations after intervention; ttest_day28, ttest results when intervention occurs on day 28 of epidemic; M_S, ttest results when matched on number of sampled susceptible individuals at time of intervention; M_SER, ttest results when matched on number of sampled non-infectious individuals at time of intervention; ttest_deg8, ttest results when average degree is 8; ttest_deg30, ttest results when average degree is 30.

Folder 3: figs (figures of paper)
- Figures created individually from code in code/generate_figures
