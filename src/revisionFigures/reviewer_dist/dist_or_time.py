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


#This figure is in response to Editor comment 1.1
#First, as you mention in lines 102-106, linkage disequilibrium results from 
#new mutations that occur in particular genetic background and then 
#disassociate with recombination. Then, what you may expect from a within-host
#population of HIV is an equilibrium level of LD due to recurrent mutation and
#recombination, particularly when one anticipates a characteristic genealogy
#(phylogeny) of within-host HIV (for example, Figure 1I of Grenfell et al. 
#2004) that shows continuous turn-over of viral lineages. I therefore expect
#that what is driving the relationship between D' ratio and distance*time in
#Figure 2 is mostly distance (at a given time). I wonder whether you actually
#see D' decrease over time while distance is fixed.


NUM_BOOTSTRAPS = 1000
# NUM_REPS = 500
NUM_GROUPS = 10

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

dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/revision/simulated_estimates_selection/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/revision/reviewer_dist/'

######################### Configure Plot Settings #############################
#plot the estimates to show how accurate they are
params = {'figure.figsize':(7.5, 3), 'axes.labelsize': 6,'axes.titlesize':6,  
          'legend.fontsize': 6, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
          'legend.title_fontsize': 6}
linewidth = 1
markerSize = 3
rcParams.update(params)
fig, axs = plt.subplots(1, 2, sharey=True)

###############################################################################
###############################################################################
###############################################################################
#Load the dataframes that will be plotted
all_stat_dfs = pd.read_pickle(dataDir + "all_stat_dfs_"+ \
                                        str(NUM_BOOTSTRAPS) + "_" + \
                                            str(NUM_GROUPS) + ".pkl")

all_stat_dfs = slim_util.label_rho(all_stat_dfs)
all_stat_dfs['Dist'] = all_stat_dfs['Locus_2'] - all_stat_dfs['Locus_1']
all_stat_time_fixed = all_stat_dfs[all_stat_dfs['Time_Diff'] == 150]
print(all_stat_dfs.columns)

############################# Run the analysis #################################
all_rho_bins = []
all_rho_bins_dist = []
all_rho_bins_time = []

#Bin the points for each rho value
for curr_rho in all_stat_dfs['Sim_Rho'].unique():
    if curr_rho not in ACC_LIST:
        continue
    curr_stat_df = all_stat_dfs[all_stat_dfs['Sim_Rho'] == curr_rho]
    curr_stat_time_fixed = all_stat_time_fixed[all_stat_time_fixed['Sim_Rho'] == curr_rho]

    #Bin the d' ratios based on distance and time
    binned_rat, binedges, bin_nums = binned_statistic(
        np.array(curr_stat_df['Dist_X_Time'].tolist()), 
        np.array(curr_stat_df['d_ratio'].tolist()), bins = 100)
    
    
    curr_rho_bins = pd.DataFrame({'Dist_X_Time': binedges[1:], 
                                  'D_Ratio': binned_rat})
    curr_rho_bins['Sim_Rho'] = curr_rho
    curr_rho_bins['Sim_float_rho'] = curr_stat_df['Sim_float_rho'].unique()[0]

    #Bin the d' ratios based on distance only
    binned_rat_dist, binedges_dist, bin_nums_dist = binned_statistic(
        np.array(curr_stat_time_fixed['Dist'].tolist()), 
        np.array(curr_stat_time_fixed['d_ratio'].tolist()), bins = 100)

    
    curr_rho_bins_dist = pd.DataFrame({'Distance': binedges_dist[1:],
                                        'D_Ratio': binned_rat_dist})
    curr_rho_bins_dist['Sim_Rho'] = curr_rho
    curr_rho_bins_dist['Sim_float_rho'] = curr_stat_df['Sim_float_rho'].unique()[0]

    #Bin the d' ratios based on time only
    # #Make sure every bin has the same number of points
    # curr_time_df  = []
    # time_counts = curr_stat_df['Time_Diff'].value_counts()
    # min_sample = min(time_counts)
    # for name, group in curr_stat_df.groupby('Time_Diff'):
    #     if len(group) > min_sample:
    #         curr_time_results = group.sample(min_sample)
    #     else:
    #         curr_time_results = group
    #     curr_time_df.append(curr_time_results)
    # curr_time_df = pd.concat(curr_time_df, ignore_index=True)
    # print(curr_time_df['Time_Diff'].value_counts())


    # binned_rat_time, binedges_time, bin_nums_time = binned_statistic(
    #     np.array(curr_stat_df['Time_Diff'].tolist()), 
    #     np.array(curr_stat_df['d_ratio'].tolist()), bins = 20)
    
    
    # curr_rho_bins_time = pd.DataFrame({'Time': binedges_time[1:], 
    #                                    'D_Ratio': binned_rat_time})
    # curr_rho_bins_time['Sim_Rho'] = curr_rho
    # curr_rho_bins_time['Sim_float_rho'] = curr_stat_df['Sim_float_rho'].unique()[0]

    all_rho_bins.append(curr_rho_bins)
    all_rho_bins_dist.append(curr_rho_bins_dist)
    # all_rho_bins_time.append(curr_rho_bins_time)

all_rho_bins = pd.concat(all_rho_bins, ignore_index=True)
all_rho_bins_dist = pd.concat(all_rho_bins_dist, ignore_index=True)
# all_rho_bins_time = pd.concat(all_rho_bins_time, ignore_index=True)

all_rho_bins['Sim_Rho'] = all_rho_bins['Sim_Rho'].astype('category')
all_rho_bins_dist['Sim_Rho'] = all_rho_bins_dist['Sim_Rho'].astype('category')
# all_rho_bins_time['Sim_Rho'] = all_rho_bins_time['Sim_Rho'].astype('category')


########################## Panel A Distance and Time ##########################
sns.lineplot(data=all_rho_bins, x='Dist_X_Time', y='D_Ratio', ax = axs[0],
             hue='Sim_Rho', linewidth = 1, 
             palette=sns.color_palette("viridis", 
             n_colors=len(all_rho_bins['Sim_Rho'].unique())), 
             hue_order = ACC_LIST)

axs[0].set_xlabel(r'Distance $\cdot$ Time (bp $\cdot$ generations)')
axs[0].set_ylabel("Linkage Decay Measure")
axs[0].legend(title = r'Simulated $\rho$')

########################## Panel B Distance ####################################
sns.lineplot(data=all_rho_bins_dist, x='Distance', y='D_Ratio', ax = axs[1],
             hue='Sim_Rho', linewidth = 1, 
             palette=sns.color_palette("viridis", 
             n_colors=len(all_rho_bins['Sim_Rho'].unique())), 
             hue_order = ACC_LIST)

axs[1].set_xlabel(r'Distance (bp)')
axs[1].set_ylabel("Linkage Decay Measure")
axs[1].get_legend().remove()

# ########################## Panel C Time ##########################
# sns.lineplot(data=all_rho_bins_time, x='Time', y='D_Ratio', ax = axs[2],
#              hue='Sim_Rho', linewidth = 1, 
#              palette=sns.color_palette("viridis", 
#              n_colors=len(all_rho_bins['Sim_Rho'].unique())), 
#              hue_order = ACC_LIST)

# axs[2].set_xlabel(r'Time (generations)')
# axs[2].set_ylabel("Linkage Decay Measure")
# axs[2].get_legend().remove()

plt.savefig(outDir + 'time_fixed_150.png', dpi = 300)


