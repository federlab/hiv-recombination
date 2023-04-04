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
NUM_BOOTSTRAPS = 1000
# NUM_BOOTSTRAPS = 10
NUM_REPS = 200
NUM_GROUPS = 100


#Today I am going to prototype an analysis to determine how well we can
#discriminate between two different recombination rates that are a fixed
#distance apart.

dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_10_03_neutral/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/croiAnalysis/'

dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_10_03_neutral/'
outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/croiAnalysis/'

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

#Filter the groups so the first d' value is greater than 0.2
all_stat_dfs = all_stat_dfs[all_stat_dfs['d_i'] > 0.2]

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
params = {'figure.figsize':(16, 6), 'axes.labelsize': 18,'axes.titlesize': 18, 
           'legend.fontsize': 18, 'xtick.labelsize': 18, 'ytick.labelsize': 18,
           'legend.title_fontsize': 18
        }
rcParams.update(params)



fig, axs = plt.subplots(1, 2)
sns.stripplot(x = 'Sim_Rho', y = 'est_rho', data = all_conf_ints, 
    jitter = True, color = 'k', s = 5, alpha = 0.3, ax = axs[1],
    order = [r"$2\times10^{-6}$", r"$10^{-5}$", r"$2\times10^{-5}$", r"$10^{-4}$", r"$2\times10^{-4}$", r"$10^{-3}$"])

# distance across the "X" or "Y" stipplot column to span, in this case 40%
label_width = 0.4

plt.tight_layout()
plt.savefig(outDir + "sim_fig_1BC_combined" + str(NUM_BOOTSTRAPS) + ".png", dpi = 300)

# params = {'figure.figsize':(3, 3), 'axes.labelsize': 8,'axes.titlesize':8,  'legend.fontsize': 8, 'xtick.labelsize': 8, 'ytick.labelsize': 8}
# rcParams.update(params)
            
for tick, text in zip(axs[1].get_xticks(), axs[1].get_xticklabels()):
    sample_name = text.get_text()  # "X" or "Y"

    #get the float value of rho corresponding with the tick
    rho_val = all_conf_ints[all_conf_ints['Sim_Rho'] == sample_name]
    rho_val = rho_val['Sim_float_rho'].unique()[0]

    # plot horizontal lines across the column, centered on the tick
    axs[1].plot([tick-label_width/2, tick+label_width/2], [rho_val, rho_val],
            lw=1, color='k')


axs[1].set_ylabel(r'Estimated Recombination Rate ($\hat{\rho}$)')
axs[1].set_xlabel(r'Simulated Recombination Rate ($\rho$)')
axs[1].set_yscale('log')


########################## Plotting Figure 1C ################################
all_rho_bins = []

#Bin the points for each rho value
for curr_rho in all_stat_dfs['Sim_Rho'].unique():
    curr_stat_df = all_stat_dfs[all_stat_dfs['Sim_Rho'] == curr_rho]

    #Bin the d' ratios so they are easier to view on the plots
    binned_rat, binedges, bin_nums = binned_statistic(
        curr_stat_df['Dist_X_Time'].to_numpy(), 
        curr_stat_df['d_ratio'].to_numpy(), bins = 100)
    
    curr_rho_bins = pd.DataFrame({'Dist_X_Time': binedges[1:], 'D_Ratio': binned_rat})
    curr_rho_bins['Sim_Rho'] = curr_rho

    all_rho_bins.append(curr_rho_bins)

all_rho_bins = pd.concat(all_rho_bins, ignore_index=True)

# #Plot our estimates against each other 
#make the rho values ints rather than strings
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
float_to_str = {0.001 : r"$10^{-3}$",
            0.0001 : r"$10^{-4}$",
            0.0002 : r"$2\times10^{-4}$",
            0.00001 : r"$10^{-5}$",
            0.00002 : r"$2\times10^{-5}$",
            0.000002 : r"$2\times10^{-6}$"}

#redo the labeling on the rho values from what was used in the simulation names
intRhoList = []
newStringRho = []
for entry in all_rho_bins['Sim_Rho']:
    intRhoList.append(rhoDict[entry])
    newStringRho.append(rho_dict_fix_strings[entry])
all_rho_bins['Sim_float_rho'] = intRhoList
all_rho_bins['Sim_Rho'] = newStringRho

sns.lineplot(data=all_rho_bins, x='Dist_X_Time', y='D_Ratio', hue='Sim_float_rho',
             linewidth = 3, ax = axs[0],
             palette=sns.color_palette("coolwarm", n_colors=len(all_rho_bins['Sim_float_rho'].unique())))
             
axs[0].legend(loc='center left', bbox_to_anchor=(1, 0.5))
axs[0].legend_.set_title(r'$\rho$')

new_labels = [r"$2\times10^{-6}$", r"$10^{-5}$", r"$2\times10^{-5}$", r"$10^{-4}$", r"$2\times10^{-4}$", r"$10^{-3}$" ]
for t, l in zip(axs[0].legend_.texts, new_labels):
    t.set_text(l)
    # t.set_font('')
axs[0].set_xlabel(r'Distance $\cdot$ Time (bp $\cdot$ generations)')
axs[0].set_ylabel("D\' Ratio")

plt.tight_layout()
plt.savefig(outDir + "sim_fig_1BC_combined" + str(NUM_BOOTSTRAPS) + ".png", dpi = 300)
plt.close()