import sys
from turtle import filling
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

#This script is going to be used to count the number of D' ratios we can expect
#from each simulation run.

dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_08_30_test/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_08_30_test/'

estimate_df = [] 
seg_hist_data = []
all_stat_df = []

#loop through each of the dataframes for the separeate simulations
for curr_data in os.listdir(dataDir):
    #only get the data directories, not hidden files
    if curr_data[0] == '.':
        continue

    #get the information for the current run
    run_info = curr_data.split('_')
    sim_rho = run_info[1]
    sim_rho = sim_rho[3:]
    rep = run_info[-1]
    sel_param = run_info[-2]

    #get the dataframe for the current run
    d_ratio_file = dataDir + curr_data + "/linkage/d_ratio"
    stat_df = pd.read_pickle(d_ratio_file)

    #label the parameters of the simulation
    stat_df['sim_rho'] = sim_rho
    stat_df['rep'] = rep
    stat_df['sel_param'] = sel_param
    
    all_stat_df.append(stat_df)

#group all of the results into one dataframe
all_stat_df = pd.concat(all_stat_df, ignore_index=True)
sim_rho_list = all_stat_df['sim_rho'].unique()
print(all_stat_df)

#count the number of segregating loci per simulation
num_seg_pairs = []
for name, group in all_stat_df.groupby(['sim_rho', 'sel_param', 'rep']):
    num_seg_pairs.append([name[0], name[1], name[2] , len(group)])

num_seg_pairs = pd.DataFrame(num_seg_pairs, columns=['sim_rho', 'sel_param', 'rep', 'num_seg_pairs'])
num_seg_pairs['sims_needed'] = np.round(np.divide(50000, num_seg_pairs['num_seg_pairs']))
print(num_seg_pairs)

# for curr_rho in sim_rho_list:
#     stat_df = all_stat_df[all_stat_df['sim_rho'] == curr_rho]

#     #get the estimate and fit for the current dataset and sample size
#     x_vals = stat_df['Dist_X_Time'].unique()
#     coeffs, fit_dat = optimize.curve_fit(plne.neher_leitner, 
#         stat_df['Dist_X_Time'].to_list(), stat_df['d_ratio'].to_list(),
#         p0 = [0, 0.26, .0000439], maxfev = 10000)
#     fit_vals = [plne.neher_leitner(x, coeffs[0], coeffs[1], coeffs[2])
#                 for x in x_vals]

#     #Bin the d' ratios so they are easier to view on the plots
#     binned_rat, binedges, bin_nums = binned_statistic(
#         stat_df['Dist_X_Time'].to_list(), 
#         stat_df['d_ratio'].to_list(), bins = 100)

#     estimate_df.append([coeffs[0], coeffs[1], coeffs[2], 
#         coeffs[1] * coeffs[2], curr_data, curr_rho, 
#         curr_data])  

#     sns.lineplot(x = binedges[:-1], y = binned_rat, color = 'tab:blue')
#     sns.lineplot(x = x_vals, y = fit_vals, color = 'tab:orange')
#     plt.xlabel(r'Distance X Time (bp/generation)')
#     plt.ylabel("D\' Ratio")


#     plt.tight_layout()
#     plt.savefig(outDir + curr_rho + 'fit.jpg')
#     plt.close()

# estimate_df = pd.DataFrame(estimate_df, columns=["C0", "C1", "C2",
#                      "Est_Rho", 'Dataset', 'Sim_Rho', 'data'] )

# num_seg_pairs = pd.DataFrame(num_seg_pairs, columns = ['Sim_Rho', 'Num_Seg_Pairs'])
# print(num_seg_pairs)