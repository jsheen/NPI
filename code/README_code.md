# Scripts used to generate simulation results

Folder 1: generate_results (scripts used to generate results of paper)
- NPI_SS_Formulae.R: script used to obtain approximate sample sizes needed.
- create_sim_bank.py (Algorithm 1): script used to create cluster simulations, with and without an enacted NPI intervention.
- create_sim_bank_day28 (Algorithm 1): script used to create cluster simulations, with and without an enacted NPI intervention on day 28 of simulation.
- search_ss_ttest.py (Algorithm 2): algorithm to find sufficient number of clusters one generation after intervention when performing a two-sample Welch's t-test (unequal variances).
- search_ss_ttest_gen2.py (Algorithm 2): algorithm to find sufficient number of clusters two generations after intervention when performing a two-sample Welch's t-test (unequal variances).
- search_ss_ttest_gen3.py (Algorithm 2): algorithm to find sufficient number of clusters three generations after intervention when performing a two-sample Welch's t-test (unequal variances).
- search_ss_M_S.py (Algorithm 2): algorithm to find sufficient number of clusters when matching on number of susceptible individuals and performing a paired one-sample t-test on the differences of the outcome.
- search_ss_M_SER.py (Algorithm 2): algorithm to find sufficient number of clusters when matching on number of non-infectious individuals and performing a paired one-sample t-test on the differences of the outcome.
- search_ss_ttest_day28 (Algorithm 2).py: algorithm to find sufficient number of clusters one generation after intervention when intervention is enacted on day 28 of the simulation when performing a two-sample Welch's t-test (unequal variances).

Folder 2: generate_tables (scripts used to generate tables of paper)
- get_ttest_rerun_csv.R: script to create csv of ttest rerun results (one generation after intervention)
- get_ttest_gen2_csv.R: script to create csv of ttest results (two generations after intervention)
- get_ttest_gen3_csv.R: script to create csv of ttest results (three generations after intervention)
- get_M_S_csv.R: script to create csv of ttest results matching on susceptible individuals
- get_M_SER_csv.R: script to create csv of ttest results matching on non-infectious individuals
- get_ttest_day28_csv.R: script to create csv of ttest results (one generation after intervention) when intervention is enacted on day 28 of the simulation.

Folder 3: generate_figures (scripts used to generate figures of paper)
- Each filename corresponds to their created figure
