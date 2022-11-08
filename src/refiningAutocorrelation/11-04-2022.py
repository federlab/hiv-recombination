import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import zaniniUtil as zu
from scipy.stats import binned_statistic
import plot_neher as plne
import autocorrelation as autocorr
from scipy import optimize
from sklearn.metrics import mean_squared_error
from matplotlib import rcParams

DIST_TIME_MAX = 50000
NUM_BOOTSTRAPS = 10
NUM_REPS = 60
NUM_GROUPS = 30

#Today I am going to prototype an analysis to determine how well we can
#discriminate between two different recombination rates that are a fixed
#distance apart.

dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_10_03_neutral/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/paper/fig2/'

#First I am going to read the data and randomly pair simulations
all_stat_dfs = []

#loop through each of the dataframes for the separate simulations
for curr_data in os.listdir(dataDir):
    #only get the data directories, not hidden files
    if curr_data[0] == '.':
        continue

    #get the information for the current run
    run_info = curr_data.split('_')
    sim_rho = run_info[1]
    sim_rho = sim_rho[3:]
    rep = run_info[-1]

    #get the dataframe for the current run
    d_ratio_file = dataDir + curr_data + "/linkage/d_ratio"
    stat_df = pd.read_pickle(d_ratio_file)
    stat_df['rep'] = int(rep[3:])
    stat_df['Sim_Rho'] = sim_rho
    all_stat_dfs.append(stat_df)
all_stat_dfs = pd.concat(all_stat_dfs)

#Randomly divide the reps into 10 groups
rep_groups = np.array(range(0, NUM_REPS))
np.random.shuffle(rep_groups)
rep_groups = np.array_split(rep_groups, NUM_GROUPS)


#Make a dictionary to label each group
group_dict = {}
for i in range(len(rep_groups)):
    for j in rep_groups[i]:
        group_dict[j] = i
        
group_labels = [group_dict[x] for x in all_stat_dfs['rep']]
all_stat_dfs['iter_group'] = group_labels


########################## Estimating recombination rates #####################
#loop through each of the distance cutoffs
all_stat_dfs = all_stat_dfs[all_stat_dfs['Dist_X_Time'] <= DIST_TIME_MAX]
all_estimate_df = []
all_fit_df = []

#loop through each rho value
for curr_rho in all_stat_dfs['Sim_Rho'].unique():
    curr_rho_stat = all_stat_dfs[all_stat_dfs['Sim_Rho'] == curr_rho]

    #loop through each iter group
    for curr_iteration in range(0, NUM_GROUPS):
        #get the data for the current rho and iteration
        curr_stat_df = curr_rho_stat[curr_rho_stat['iter_group'] == curr_iteration]

        #Get the current estimate
        lower_fit, upper_fit, estimate_df = plne.bootstrap_rho(curr_stat_df,
                                                             NUM_BOOTSTRAPS)
        estimate_df['Group'] = curr_iteration
        estimate_df['Sim_Rho'] = curr_rho
        all_estimate_df.append(estimate_df)


all_estimate_df = pd.concat(all_estimate_df, ignore_index=True)

all_conf_ints = []

#Calculate the confidence intervals for the estimates with low + high VL
for name, group in all_estimate_df.groupby(['Group', 'Sim_Rho']):
    lower_conf = np.quantile(group['Estimated_Rho'], 0.025)
    mid_est = np.quantile(group['Estimated_Rho'], 0.5)
    upper_conf = np.quantile(group['Estimated_Rho'], 0.975)

    #append the confidence interval
    all_conf_ints.append([lower_conf, mid_est, upper_conf, name[0], name[1]])

all_conf_ints = pd.DataFrame(all_conf_ints, 
    columns=['lower_conf', 'est_rho', 'upper_conf', 'Group', 'Sim_Rho'])

print(all_conf_ints)

########################## Comparing pairs of groups ##########################
#To get the comparison, we need to convert the rho values to floats
rhoDict = {"0.001" : 0.001,
            "1e-04" : 0.0001,
            "2e-04" : 0.0002,
            "1e-05" : 0.00001,
            "2e-05" : 0.00002,
            "2e-06" : 0.000002}
rho_dict_fix_strings = { "0.001" : r"$10^{-3}$",
                        "1e-04" : r"$10^{-4}$",
                        "2e-04" : r"$2\times10^{-4}$",
                        "1e-05" : r"$10^{-5}$",
                        "2e-05" : r"$2\times10^{-5}$",
                        "2e-06" : r"$2\times10^{-6}$"}

#redo the labeling on the rho values from what was used in the simulation names
intRhoList = []
newStringRho = []
for entry in all_conf_ints['Sim_Rho']:
    intRhoList.append(rhoDict[entry])
    newStringRho.append(rho_dict_fix_strings[entry])
all_conf_ints['Sim_float_rho'] = intRhoList
all_conf_ints['Sim_Rho'] = newStringRho


#We need to loop through all pairs of rates
possible_rhos = all_conf_ints['Sim_float_rho'].unique()
possible_rhos = np.sort(possible_rhos)

pair_results_df = []

#loop through all but the highest rho value
for i in range(len(possible_rhos -1)):
    curr_rho_1 = possible_rhos[i]
    rho_1_df = all_conf_ints[all_conf_ints['Sim_float_rho'] == curr_rho_1]
    curr_str_rho_1 = rho_1_df['Sim_Rho'].unique()[0]

    for j in range(i+1, len(possible_rhos)):
        curr_rho_2 = possible_rhos[j]
        rho_2_df = all_conf_ints[all_conf_ints['Sim_float_rho'] == curr_rho_2]
        curr_str_rho_2 = rho_2_df['Sim_Rho'].unique()[0]

        pair_diff = curr_rho_2/curr_rho_1

        #randomly pair the simulations
        rep_list_1 = rho_1_df['Group'].unique()
        np.random.shuffle(rep_list_1)
        rep_list_2 = rho_2_df['Group'].unique()
        np.random.shuffle(rep_list_2)
        assert len(rep_list_1) == len(rep_list_2), "The number of replicates for each rho value is not the same"
        rep_pairs = [(rep_list_1[i], rep_list_2[i]) for i in range(len(rep_list_1))]
        
        #loop through each pair and check if the confidence intervals overlap
        for curr_pair in rep_pairs:
            rho_1_group = rho_1_df[rho_1_df['Group'] == curr_pair[0]]
            rho_2_group = rho_2_df[rho_2_df['Group'] == curr_pair[1]]

            #check if the estimates are in the correct order
            if rho_1_group['est_rho'].values[0] > rho_2_group['est_rho'].values[0]:
                order_correct = False
            else: order_correct = True  
            
            #check if the confidence intervals overlap
            lower_conf_1 = rho_1_group['lower_conf'].values[0]
            upper_conf_1 = rho_1_group['upper_conf'].values[0]
            lower_conf_2 = rho_2_group['lower_conf'].values[0]
            upper_conf_2 = rho_2_group['upper_conf'].values[0]
            
            if (lower_conf_1 < lower_conf_2) and (lower_conf_2 < upper_conf_1):
                overlap = True
            elif (lower_conf_1 < upper_conf_2) and (upper_conf_2 < upper_conf_1):
                overlap = True
            else:
                overlap = False

            #append the results
            pair_results_df.append([curr_rho_1, curr_rho_2, order_correct, overlap, pair_diff])

pair_results_df = pd.DataFrame(pair_results_df, columns=['rho_1', 'rho_2', 'order_correct', 'overlap', 'pair_diff'])

#Now calculate the proportion correct 

########################## Plotting our Results ###############################          

sns.stripplot(x = 'Sim_Rho', y = 'Est_Rho', data = estimate_df, 
    jitter = True, color = 'k', s = 8, alpha = 0.3,
    order = [r"$2\times10^{-6}$", r"$10^{-5}$", r"$2\times10^{-5}$", r"$10^{-4}$", r"$2\times10^{-4}$", r"$10^{-3}$"])