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
from scipy.stats import binned_statistic


DIST_TIME_MAX = 50000
RATIOS_PER_GROUP = 25000

NUM_BOOTSTRAPS = 1000
NUM_REPS = 200
NUM_GROUPS = 40


#Today I am reworking the simulation figure to use formatting that is more
#clear to people.

dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/paper/simulated_estimates/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_10_03_neutral/04-20-2023/'

# dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/paper/simulated_estimates/'
# outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_10_03_neutral/04-20-2023/'
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
all_stat_dfs = pd.read_pickle(dataDir + "all_stat_dfs_"+ \
                                        str(NUM_BOOTSTRAPS) + "_" + \
                                            str(NUM_GROUPS) + ".pkl")

all_stat_dfs = label_rho(all_stat_dfs)

########################## Plotting Figure 1C ################################
params = {'figure.figsize':(4.5, 2.5), 'axes.labelsize': 8,'axes.titlesize':8,  'legend.fontsize': 8, 'xtick.labelsize': 8, 'ytick.labelsize': 8,
          'legend.title_fontsize': 8}
rcParams.update(params)

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
    curr_rho_bins['Sim_float_rho'] = curr_stat_df['Sim_float_rho'].unique()[0]

    all_rho_bins.append(curr_rho_bins)

all_rho_bins = pd.concat(all_rho_bins)

fig = sns.lineplot(data=all_rho_bins, x='Dist_X_Time', y='D_Ratio', hue='Sim_float_rho',
             linewidth = 1,
             palette=sns.color_palette("viridis", n_colors=len(all_rho_bins['Sim_float_rho'].unique())))

fig.legend(loc='center left', bbox_to_anchor=(1, 0.5),title = r'Recombination Rate ($\rho$)')
# fig.legend_.set_title(r'$\rho$')


new_labels = [r"$2\times10^{-6}$", r"$10^{-5}$", r"$2\times10^{-5}$", r"$10^{-4}$", r"$2\times10^{-4}$", r"$10^{-3}$" ]
for t, l in zip(fig.legend_.texts, new_labels):
    t.set_text(l)
    # t.set_font('')
plt.xlabel(r'Distance $\cdot$ Time (bp $\cdot$ generations)')
plt.ylabel("D\' Ratio")
plt.tight_layout()
plt.savefig(outDir + "fig2_option2.png", dpi = 300)
plt.close()