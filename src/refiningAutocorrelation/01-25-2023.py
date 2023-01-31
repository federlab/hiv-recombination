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
from matplotlib import rcParams

DIST_TIME_MAX = 50000
NUM_BOOTSTRAPS = 1000
NUM_BOOTSTRAPS = 10
NUM_REPS = 20#0
NUM_GROUPS = 2
DI_THRESH = 0.2


#Today I am going to try to look at the distribution of D_i vs distance and time

dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_10_03_neutral/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_10_03_neutral/01-25-2023/'

# dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_10_03_neutral/'
# outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_10_03_neutral/01-25-2023/'

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
all_stat_dfs = all_stat_dfs[all_stat_dfs['d_i'] > DI_THRESH]
print(all_stat_dfs.head())

#Filter the groups so the first d' value is greater than 0.2
# all_stat_dfs = all_stat_dfs[all_stat_dfs['d_i'] > 0.2]
all_stat_dfs = all_stat_dfs[all_stat_dfs['Dist_X_Time'] < DIST_TIME_MAX]

my_col_order = ['0.001', "1e-04", "2e-04", "1e-05", "2e-05", "2e-06"]

#Print a scatter plot of initial linkage vs distance and time
my_grid = sns.FacetGrid(all_stat_dfs, col ='Sim_Rho', col_wrap = 3, col_order = my_col_order )
print(all_stat_dfs['Dist_X_Time'])
my_grid.map(sns.scatterplot, 'Dist_X_Time', 'd_i', alpha = 0.0075)
plt.savefig(outDir + 'd_i_vs_distTime.png')

#Print a scatter plot of initial linkage vs distance and time
my_grid = sns.FacetGrid(all_stat_dfs, col ='Sim_Rho', col_wrap = 3, col_order = my_col_order )
print(all_stat_dfs['Dist_X_Time'])
my_grid.map(sns.histplot, 'Dist_X_Time')
plt.ylim(0, 15000)
plt.savefig(outDir + 'distTime.png')

#Print a histogram of initial linkage values
my_grid = sns.FacetGrid(all_stat_dfs, col ='Sim_Rho', col_wrap = 3, col_order = my_col_order )
print(all_stat_dfs['Dist_X_Time'])
my_grid.map(sns.histplot, 'd_i')
plt.savefig(outDir + 'd_i_hist.png')

#Now repeat the process with a threshold of 0.2
all_stat_dfs = all_stat_dfs[all_stat_dfs['d_i'] > DI_THRESH]

#Print a scatter plot of initial linkage vs distance and time
my_grid = sns.FacetGrid(all_stat_dfs, col ='Sim_Rho', col_wrap = 3, col_order = my_col_order )
print(all_stat_dfs['Dist_X_Time'])
my_grid.map(sns.histplot, 'Dist_X_Time')
plt.ylim(0, 15000)
plt.savefig(outDir + 'distTime_'+ str(DI_THRESH) + '_.png')


#Now I am going to try to downsample the segregating loci by bin
#Make a dictionary to label each group
group_dict = {}
for i in range(len(rep_groups)):
    for j in rep_groups[i]:
        group_dict[j] = i
        
group_labels = [group_dict[x] for x in all_stat_dfs['rep']]
all_stat_dfs['iter_group'] = group_labels


#loop through each of the sample sizes
estimate_df_size = []
for curr_size in D_PRIME_NUMS:
    print("Size is: " + str(curr_size))
    #loop through each rho value
    for curr_rho in all_stat_dfs['Sim_Rho'].unique():
        curr_rho_stat = all_stat_dfs[all_stat_dfs['Sim_Rho'] == curr_rho]
        
        #Add the distance example cutoff
        curr_rho_stat = curr_rho_stat[curr_rho_stat['Dist_X_Time'] <= DIST_EXAMPLE] 

        #loop through each iter group
        for curr_iteration in range(1, NUM_GROUPS+1):
            #get the data for the current rho and iteration
            curr_stat_df = curr_rho_stat[curr_rho_stat['rep'] == curr_iteration]

            #make 300 bins

            #sample segregating sites while the D' num is below the threshold
            sample_df = []
            sample_size = 0
            sampled_loci = set()
            unsampled_loci = set(list(curr_stat_df['Locus_1'].unique()) + list(curr_stat_df['Locus_2'].unique().tolist()))
            while sample_size < curr_size and len(unsampled_loci) > 0:
                #sample another locus
                my_sample = np.random.choice(list(unsampled_loci))
                sampled_loci.add(my_sample)
                unsampled_loci.remove(my_sample)
                sample_df = curr_stat_df[curr_stat_df['Locus_1'].isin(sampled_loci) & (curr_stat_df['Locus_2'].isin(sampled_loci))]
                sample_size = len(sample_df)
                