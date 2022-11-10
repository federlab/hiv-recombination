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

THRESHOLD = 0.2
DIST_TIME_MAX = 50000
NUM_REPS = 60

#I am going to try pooling all of the neutral data and plotting it to show how the decay depends on recombination rate

dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_10_03_neutral/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_10_03_neutral/'

#First, we will get all of the data and divide it into groups
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
    if int(rep[3:]) > NUM_REPS:
        continue

    #get the dataframe for the current run
    d_ratio_file = dataDir + curr_data + "/linkage/d_ratio"
    stat_df = pd.read_pickle(d_ratio_file)
    stat_df['rep'] = int(rep[3:])
    stat_df['Sim_Rho'] = sim_rho
    all_stat_dfs.append(stat_df)
all_stat_dfs = pd.concat(all_stat_dfs)

all_stat_dfs = all_stat_dfs[all_stat_dfs['Dist_X_Time'] <= DIST_TIME_MAX]

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
sns.lineplot(data=all_rho_bins, x='Dist_X_Time', y='D_Ratio', hue='Sim_Rho',
             palette=sns.color_palette("icefire", n_colors=len(all_stat_dfs['Sim_Rho'].unique())))
plt.savefig(outDir + 'rho_dist_x_time.png')