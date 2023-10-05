import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import slimUtil as slim_util
from scipy.stats import binned_statistic
from matplotlib import rcParams
import plot_neher as plne
from scipy.stats import linregress

#Make a 9 panel figure of the selection curves fits faceted by rho
dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/revision/simulated_estimates_selection/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/revision/supp_figs/supp_sel_curves/'

ACC_LIST = [r"$2\times10^{-6}$",
                        r"$5\times10^{-6}$",
                        r"$10^{-5}$",
                        r"$2\times10^{-5}$",
                        r"$5\times10^{-5}$",
                        r"$10^{-4}$",
                        r"$2\times10^{-4}$",
                        r"$5\times10^{-4}$",
                        r"$10^{-3}$"]
NUM_BOOTSTRAPS = 1000
NUM_GROUPS = 10


######################### Set the plotting parameters #########################
rcParams.update({'figure.figsize':(7,8), 'axes.labelsize': 6,
                'axes.titlesize':6, 'legend.fontsize': 6, 
                'xtick.labelsize': 6, 'ytick.labelsize': 6,
                'legend.title_fontsize': 6})
fig, axs = plt.subplots(3, 3, sharex = True, sharey = True)
iter_ax = axs.flatten()

###############################################################################

#Load the dataframes that will be plotted
all_stat_dfs = pd.read_pickle(dataDir + "all_stat_dfs_"+ \
                                        str(NUM_BOOTSTRAPS) + "_" + \
                                            str(NUM_GROUPS) + ".pkl")

all_stat_dfs = slim_util.label_rho(all_stat_dfs)
all_stat_dfs = all_stat_dfs[all_stat_dfs['Sim_Rho'].isin(ACC_LIST)]

#Choose a group to plot and initialize the plot
example_group = np.random.choice(range(0, NUM_GROUPS))

all_stat_dfs = all_stat_dfs[all_stat_dfs['iter_group'] == example_group]

for i in range(len(ACC_LIST)):
    curr_rho = ACC_LIST[i]
    curr_ax = iter_ax[i]
    curr_stat_df = all_stat_dfs[all_stat_dfs['Sim_Rho'] == curr_rho]

    #Bin the d' ratios so they are easier to view on the plots
    binned_rat, binedges, bin_nums = binned_statistic(
        curr_stat_df['Dist_X_Time'].to_numpy(), 
        curr_stat_df['d_ratio'].to_numpy(), bins = 100)

    sns.lineplot(y = binned_rat, x = binedges[1:], color = 'black', 
                 ax = curr_ax)
    curr_ax.set_title(r'$\rho$ = ' + str(curr_rho), fontsize = 8)
    curr_ax.set_ylabel('Linkage Decay Measure')
    curr_ax.set_xlabel(r'Distance $\cdot$ Time (bp $\cdot$ generations)')

plt.savefig(outDir + 'sel_curves.jpg', dpi = 300)
print(all_stat_dfs.head())