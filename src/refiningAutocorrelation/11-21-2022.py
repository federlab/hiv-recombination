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
NUM_GROUPS = 10

#Today I am going to prototype an analysis to determine how well we can
#discriminate between two different recombination rates that are a fixed
#distance apart.

dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_09_07_MPL/'

# dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/slimDatasets/'
# outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_09_07_MPL/'
dataDirList = ['2022_09_07_MPL_1e-3/', '2022_09_07_MPL_1e-4/'] #, '2022_09_07_MPL_1e-5/']
dataDirReps = [160, 3270] #, 5490]

#First I am going to read the data and randomly pair simulations
all_stat_dfs = []


#loop through each of the dataframes for the separate simulations
for i in range(len(dataDirList)):
    curr_dir = dataDir + dataDirList[i]
    curr_reps = dataDirReps[i]
    rho_stat_dfs = []
    for curr_data in os.listdir(curr_dir):
        #only get the data directories, not hidden files
        if curr_data[0] == '.':
            continue

        #get the information for the current run
        run_info = curr_data.split('_')
        sim_rho = run_info[1]
        sim_rho = sim_rho[3:]
        rep = run_info[-1]

        #get the dataframe for the current run
        d_ratio_file = curr_dir + curr_data + "/linkage/d_ratio"
        stat_df = pd.read_pickle(d_ratio_file)
        stat_df['rep'] = int(rep[3:])
        stat_df['Sim_Rho'] = sim_rho
        rho_stat_dfs.append(stat_df)
    rho_stat_dfs = pd.concat(rho_stat_dfs)

    #Randomly divide the reps into 10 groups
    rep_groups = np.array(range(0, curr_reps))
    np.random.shuffle(rep_groups)
    rep_groups = np.array_split(rep_groups, NUM_GROUPS)


    #Make a dictionary to label each group
    group_dict = {}
    for i in range(len(rep_groups)):
        for j in rep_groups[i]:
            group_dict[j] = i
            
    group_labels = [group_dict[x] for x in rho_stat_dfs['rep']]
    rho_stat_dfs['iter_group'] = group_labels
    all_stat_dfs.append(rho_stat_dfs)

all_stat_dfs = pd.concat(all_stat_dfs, ignore_index=True)


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

########################## Plotting the accuracy of the estimates #############
#plot the estimates to show how accurate they are
sns.set(rc={'figure.figsize':(10,10)}, font_scale = 2, font = '')
sns.set_palette("tab10")
sns.set_style("white")

fig = sns.stripplot(x = 'Sim_Rho', y = 'est_rho', data = all_conf_ints, 
    jitter = True, color = 'k', s = 8, alpha = 0.3,
    order = [r"$10^{-4}$", r"$10^{-3}$"])
    #order = [r"$2\times10^{-6}$", r"$10^{-5}$", r"$2\times10^{-5}$", r"$10^{-4}$", r"$2\times10^{-4}$", r"$10^{-3}$"])

# distance across the "X" or "Y" stipplot column to span, in this case 40%
label_width = 0.4

#need to save the fig to force the axis to be drawn
plt.savefig(outDir + "selection_accuracy.png", dpi = 300)
            
for tick, text in zip(fig.get_xticks(), fig.get_xticklabels()):
    sample_name = text.get_text()  # "X" or "Y"

    #get the float value of rho corresponding with the tick
    rho_val = all_conf_ints[all_conf_ints['Sim_Rho'] == sample_name]
    rho_val = rho_val['Sim_float_rho'].unique()[0]

    # plot horizontal lines across the column, centered on the tick
    fig.plot([tick-label_width/2, tick+label_width/2], [rho_val, rho_val],
            lw=2, color='k')

plt.ylabel(r'Estimated Value of $\rho$')
plt.xlabel(r'Simulation Value of $\rho$')
plt.yscale('log')
plt.savefig(outDir + "selection_accuracy.png", dpi = 300)
plt.close()

####################### Plotting the Ddelta t distribution ####################
myfacets = sns.FacetGrid(all_stat_dfs, col='Sim_Rho', col_wrap=3, height=5, aspect=1.5)
myfacets.map(sns.histplot, 'Dist_X_Time', kde = True,  stat = 'density', bins = 300)
plt.savefig(outDir + "Ddelta_t_dist.png", dpi = 300)