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
NUM_REPS = 200
NUM_GROUPS = 10
DI_THRESHOLD_LIST = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

#Today I am going to test the accuracy  when I move up the initial linkage threshold
dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_10_03_neutral/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_10_03_neutral/'

# dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_10_03_neutral/'
# outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_10_03_neutral/'

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
    if int(rep[3:]) >= NUM_REPS:
        continue

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


print(len(all_stat_dfs))
########################## Estimating recombination rates #####################
#loop through each of the distance cutoffs
all_stat_dfs = all_stat_dfs[all_stat_dfs['Dist_X_Time'] <= DIST_TIME_MAX]
all_estimate_df = []
all_fit_df = []

#loop through each threshold
for curr_threshold in DI_THRESHOLD_LIST:
    curr_thresh_stat = all_stat_dfs[all_stat_dfs['d_i'] >= curr_threshold]
    #loop through each rho value
    for curr_rho in curr_thresh_stat['Sim_Rho'].unique():
        curr_rho_stat = curr_thresh_stat[curr_thresh_stat['Sim_Rho'] == curr_rho]

        #loop through each iter group
        for curr_iteration in range(0, NUM_GROUPS):
            #get the data for the current rho and iteration
            curr_stat_df = curr_rho_stat[curr_rho_stat['iter_group'] == curr_iteration]

            #Get the current estimate
            lower_fit, upper_fit, estimate_df = plne.bootstrap_rho(curr_stat_df,
                                                                NUM_BOOTSTRAPS)
            estimate_df['Group'] = curr_iteration
            estimate_df['Sim_Rho'] = curr_rho
            estimate_df['Threshold'] = curr_threshold
            all_estimate_df.append(estimate_df)


all_estimate_df = pd.concat(all_estimate_df, ignore_index=True)

all_conf_ints = []

#Calculate the confidence intervals for the estimates with low + high VL
for name, group in all_estimate_df.groupby(['Group', 'Sim_Rho', 'Threshold']):
    lower_conf = np.quantile(group['Estimated_Rho'], 0.025)
    mid_est = np.quantile(group['Estimated_Rho'], 0.5)
    upper_conf = np.quantile(group['Estimated_Rho'], 0.975)

    #append the confidence interval
    all_conf_ints.append([lower_conf, mid_est, upper_conf, name[0], name[1], name[2]])

all_conf_ints = pd.DataFrame(all_conf_ints, 
    columns=['lower_conf', 'est_rho', 'upper_conf', 'Group', 'Sim_Rho', 'Threshold'])


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
sns.set(rc={'figure.figsize':(35,10)}, font_scale = 2, font = '')
sns.set_palette("tab10")
sns.set_style("white")

fig, axs = plt.subplots(2, 4, sharex = True, sharey = True)
axs_nums = [(0,0),(0,1), (0,2), (0,3), (1,0), (1,1), (1,2), (1,3)]

for i in range(len(all_conf_ints['Threshold'].unique())):
    curr_plot_conf = all_conf_ints[all_conf_ints['Threshold'] == all_conf_ints['Threshold'].unique()[i]]
    print(curr_plot_conf)
    ax_coords = axs_nums[i]
    curr_ax = axs[ax_coords[0]][ax_coords[1]]

    sns.stripplot(x = 'Sim_Rho', y = 'est_rho', data = curr_plot_conf, ax = curr_ax,
        jitter = True, color = 'k', s = 8, alpha = 0.3,
        order = [r"$2\times10^{-6}$", r"$10^{-5}$", r"$2\times10^{-5}$", r"$10^{-4}$", r"$2\times10^{-4}$", r"$10^{-3}$"])
    plt.savefig(outDir + "accuracy_facet_threshold.png")
    # distance across the "X" or "Y" stipplot column to span, in this case 40%
    label_width = 0.4
                
    for tick, text in zip(curr_ax.get_xticks(), curr_ax.get_xticklabels()):
        sample_name = text.get_text()  # "X" or "Y"
        print(sample_name)

        #get the float value of rho corresponding with the tick
        rho_val = curr_plot_conf[curr_plot_conf['Sim_Rho'] == sample_name]
        rho_val = rho_val['Sim_float_rho'].unique()[0]

        # plot horizontal lines across the column, centered on the tick
        curr_ax.plot([tick-label_width/2, tick+label_width/2], [rho_val, rho_val],
                lw=2, color='k')
        axs[0, ax_coords[1]].plot([tick-label_width/2, tick+label_width/2], [rho_val, rho_val],
                lw=2, color='k')


    curr_ax.set_title("Threshold = " + str(all_conf_ints['Threshold'].unique()[i]))

axs[0][0].set_ylabel(r'Estimated Value of $\rho$')
axs[1][0].set_ylabel(r'Estimated Value of $\rho$')
axs[1][0].set_xlabel(r'Simulation Value of $\rho$')
axs[1][1].set_xlabel(r'Simulation Value of $\rho$')
axs[1][2].set_xlabel(r'Simulation Value of $\rho$')
axs[1][3].set_xlabel(r'Simulation Value of $\rho$')
plt.savefig(outDir + "accuracy_facet_threshold_unlogged.png", dpi = 300)
plt.tight_layout()
plt.yscale('log')
#plt.savefig(outDir + "accuracy_facet_threshold.png", dpi = 300)
plt.close()
