import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import estimation_util as est_util
import slimUtil as slim_util
from matplotlib import rcParams
from matplotlib.lines import Line2D
from scipy.stats import binned_statistic


#This figure is in response to Rm comment 1.3
#I am just plotting 3A as a function of \rho d \Delta t

NUM_BOOTSTRAPS = 1000
NUM_REPS = 100
NUM_GROUPS = 20

#The rho values to show in the accuracy panel
ACC_LIST = [r"$2\times10^{-6}$",
                        r"$5\times10^{-6}$",
                        r"$10^{-5}$",
                        r"$2\times10^{-5}$",
                        r"$5\times10^{-5}$",
                        r"$10^{-4}$",
                        r"$2\times10^{-4}$",
                        r"$5\times10^{-4}$",
                        r"$10^{-3}$"]

dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/revision/simulated_estimates_neutral/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/revision/supp_figs/supp_collapse/'

######################### Configure Plot Settings #############################
#plot the estimates to show how accurate they are
params = {'figure.figsize':(3, 3), 'axes.labelsize': 6,'axes.titlesize':6,  
          'legend.fontsize': 6, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
          'legend.title_fontsize': 6}
linewidth = 1
markerSize = 3
rcParams.update(params)

###############################################################################
###############################################################################
###############################################################################
#Load the dataframes that will be plotted
all_conf_ints = pd.read_pickle(dataDir + "all_conf_ints_"+ \
                                        str(NUM_BOOTSTRAPS) + "_" + \
                                            str(NUM_GROUPS) + ".pkl")
all_stat_dfs = pd.read_pickle(dataDir + "all_stat_dfs_"+ \
                                        str(NUM_BOOTSTRAPS) + "_" + \
                                            str(NUM_GROUPS) + ".pkl")

all_conf_ints = slim_util.label_rho(all_conf_ints)
all_stat_dfs = slim_util.label_rho(all_stat_dfs)

######################## Bin the values for the D' Ratio Curves ###############
all_rho_bins = []

#Bin the points for each rho value
for curr_rho in all_stat_dfs['Sim_Rho'].unique():
    if curr_rho not in ACC_LIST:
        continue
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

all_rho_bins['Sim_Rho'] = all_rho_bins['Sim_Rho'].astype('category')
all_rho_bins['RDT'] = all_rho_bins['Dist_X_Time'] * all_rho_bins['Sim_float_rho']

########################## Panel A  D' Ratio ###################################
sns.lineplot(data=all_rho_bins, x='RDT', y='D_Ratio', 
             hue='Sim_Rho', linewidth = 1, alpha = 1,
             palette=sns.color_palette("viridis", 
             n_colors=len(all_rho_bins['Sim_Rho'].unique())), 
             hue_order = ACC_LIST)

plt.xlabel(r'Recombination Rate ($\rho$) $\cdot$ Distance $\cdot$ Time')
plt.ylabel("Linkage Decay Measure")
plt.legend(title = r'Simulated $\rho$')
plt.savefig(outDir + 'collapsed_curves.png', dpi = 300)
