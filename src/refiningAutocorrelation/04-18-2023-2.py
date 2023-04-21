import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import estimation_util as est_util
from matplotlib import rcParams
from matplotlib.lines import Line2D


DIST_TIME_MAX = 50000
RATIOS_PER_GROUP = 25000

NUM_BOOTSTRAPS = 1000
NUM_REPS = 200
NUM_GROUPS = 40


#Today I am reworking the simulation figure to use formatting that is more
#clear to people.

dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/paper/simulated_estimates/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_10_03_neutral/04-18-2023/'

# dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/paper/simulated_estimates/'
# outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_10_03_neutral/04-18-2023/'
###############################################################################
#################################### Helper Functions #########################
###############################################################################

def label_rho(my_dataframe):
    #This function takes a dataframe and labels the rho values with floats
    #instead of strings so that they can be used as a continuous variable
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
    for entry in my_dataframe['Sim_Rho']:
        intRhoList.append(rhoDict[entry])
        newStringRho.append(rho_dict_fix_strings[entry])
    my_dataframe['Sim_float_rho'] = intRhoList
    my_dataframe['Sim_Rho'] = newStringRho

    return my_dataframe
###############################################################################
###############################################################################
###############################################################################

all_conf_ints = pd.read_pickle(dataDir + "all_conf_ints_"+ \
                                        str(NUM_BOOTSTRAPS) + "_" + \
                                            str(NUM_GROUPS) + ".pkl")

all_conf_ints = label_rho(all_conf_ints)



########################## Plotting the accuracy of the estimates #############
#plot the estimates to show how accurate they are
params = {'figure.figsize':(6.5, 3), 'axes.labelsize': 8,'axes.titlesize':8,  'legend.fontsize': 8, 'xtick.labelsize': 8, 'ytick.labelsize': 8}
rcParams.update(params)



fig, axs = plt.subplots(1, 2)
sns.stripplot(x = 'Sim_Rho', y = 'est_rho', data = all_conf_ints, 
    jitter = True, color = 'k', s = 3, alpha = 0.3, ax = axs[1],
    order = [r"$2\times10^{-6}$", r"$10^{-5}$", r"$2\times10^{-5}$", r"$10^{-4}$", r"$2\times10^{-4}$", r"$10^{-3}$"])

# distance across the "X" or "Y" stipplot column to span, in this case 40%
label_width = 0.5

# we need to save the figure once here so that we can access the axis labels
plt.tight_layout()
plt.savefig(outDir + "sim_fig_1BC_combined" + str(NUM_BOOTSTRAPS) + ".png", dpi = 300)

params = {'figure.figsize':(3, 3), 'axes.labelsize': 8,'axes.titlesize':8,  'legend.fontsize': 8, 'xtick.labelsize': 8, 'ytick.labelsize': 8}
rcParams.update(params)
            
for tick, text in zip(axs[1].get_xticks(), axs[1].get_xticklabels()):
    sample_name = text.get_text()  # "X" or "Y"

    #get the float value of rho corresponding with the tick
    rho_val = all_conf_ints[all_conf_ints['Sim_Rho'] == sample_name]
    rho_val = rho_val['Sim_float_rho'].unique()[0]
    estimate = np.mean(all_conf_ints[all_conf_ints['Sim_Rho'] == sample_name]['est_rho'])

    # plot horizontal lines across the column, centered on the tick
    axs[1].plot([tick-label_width/2, tick+label_width/2], [rho_val, rho_val],
            lw=1, color='r', linestyle = '--')
    axs[1].plot([tick-label_width/2, tick+label_width/2], [estimate, estimate],
            lw=1, color='k')


# make a legend showing which line is the mean and which line is the truth
redline = Line2D([0], [0], color='r', linestyle = '--', label = r'True $\rho$')
blackline = Line2D([0], [0], color='k', label = r'Mean Estimate')
axs[1].legend(handles = [redline, blackline], loc = 'upper left', frameon = True)


axs[1].set_ylabel(r'Estimated Recombination Rate ($\hat{\rho}$)')
axs[1].set_xlabel(r'Simulated Recombination Rate ($\rho$)')
axs[1].set_yscale('log')


# ########################## Plotting Figure 1C ################################
# all_rho_bins = []

# #Bin the points for each rho value
# for curr_rho in all_stat_dfs['Sim_Rho'].unique():
#     curr_stat_df = all_stat_dfs[all_stat_dfs['Sim_Rho'] == curr_rho]

#     #Bin the d' ratios so they are easier to view on the plots
#     binned_rat, binedges, bin_nums = binned_statistic(
#         curr_stat_df['Dist_X_Time'].to_numpy(), 
#         curr_stat_df['d_ratio'].to_numpy(), bins = 100)
    
#     curr_rho_bins = pd.DataFrame({'Dist_X_Time': binedges[1:], 'D_Ratio': binned_rat})
#     curr_rho_bins['Sim_Rho'] = curr_rho

#     all_rho_bins.append(curr_rho_bins)

# all_rho_bins = pd.concat(all_rho_bins, ignore_index=True)

# # #Plot our estimates against each other 
# #make the rho values ints rather than strings
# rhoDict = {"0.001" : 0.001,
#             "1e-04" : 0.0001,
#             "2e-04" : 0.0002,
#             "1e-05" : 0.00001,
#             "2e-05" : 0.00002,
#             "2e-06" : 0.000002}
# rho_dict_fix_strings = { "0.001" : r"$10^{-3}$",
#                         "1e-04" : r"$10^{-4}$",
#                         "2e-04" : r"$2\times10^{-4}$",
#                         "1e-05" : r"$10^{-5}$",
#                         "2e-05" : r"$2\times10^{-5}$",
#                         "2e-06" : r"$2\times10^{-6}$"}
# float_to_str = {0.001 : r"$10^{-3}$",
#             0.0001 : r"$10^{-4}$",
#             0.0002 : r"$2\times10^{-4}$",
#             0.00001 : r"$10^{-5}$",
#             0.00002 : r"$2\times10^{-5}$",
#             0.000002 : r"$2\times10^{-6}$"}

# #redo the labeling on the rho values from what was used in the simulation names
# intRhoList = []
# newStringRho = []
# for entry in all_rho_bins['Sim_Rho']:
#     intRhoList.append(rhoDict[entry])
#     newStringRho.append(rho_dict_fix_strings[entry])
# all_rho_bins['Sim_float_rho'] = intRhoList
# all_rho_bins['Sim_Rho'] = newStringRho

# sns.lineplot(data=all_rho_bins, x='Dist_X_Time', y='D_Ratio', hue='Sim_float_rho',
#              linewidth = 1, ax = axs[0],
#              palette=sns.color_palette("coolwarm", n_colors=len(all_rho_bins['Sim_float_rho'].unique())))

# axs[0].legend_.set_title(r'$\rho$')

# new_labels = [r"$2\times10^{-6}$", r"$10^{-5}$", r"$2\times10^{-5}$", r"$10^{-4}$", r"$2\times10^{-4}$", r"$10^{-3}$" ]
# for t, l in zip(axs[0].legend_.texts, new_labels):
#     t.set_text(l)
#     # t.set_font('')
# axs[0].set_xlabel(r'Distance $\cdot$ Time (bp $\cdot$ generations)')
# axs[0].set_ylabel("D\' Ratio")
plt.tight_layout()
plt.savefig(outDir + "sim_fig_1BC_combined" + str(NUM_BOOTSTRAPS) + "reps_"+ str(NUM_REPS) +"groups_" + str(NUM_GROUPS) + ".png", dpi = 300)
plt.close()