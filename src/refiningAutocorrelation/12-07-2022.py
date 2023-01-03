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
NUM_GROUPS = 10
NUM_REPS = 20
 

#Today I am going to plot the neutral recombination rate curves next to the selection recombination rate curves

dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/'
neutralDataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_10_03_neutral/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_09_07_MPL/'

dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/slimDatasets/'
neutralDataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_10_03_neutral/'
outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_09_07_MPL/'

dataDirList = ['2022_09_07_MPL_1e-3/', '2022_09_07_MPL_1e-4/', '2022_09_07_MPL_1e-5/']
dataDirReps = [160, 3270, 5490]

#First I am going to read the data and randomly pair simulations
all_stat_dfs = []


#loop through each of the dataframes for the separate simulations
for i in range(len(dataDirList)):
    curr_dir = dataDir + dataDirList[i]
    curr_reps = dataDirReps[i]
    rho_stat_dfs = []
    for curr_data in os.listdir(curr_dir):
        #only get the data directories, not hidden files
        if curr_data[0] == '.':
            continue

        #get the information for the current run
        run_info = curr_data.split('_')
        sim_rho = run_info[1]
        sim_rho = sim_rho[3:]
        rep = run_info[-1]

        #get the dataframe for the current run
        d_ratio_file = curr_dir + curr_data + "/linkage/d_ratio"
        stat_df = pd.read_pickle(d_ratio_file)
        stat_df['rep'] = int(rep[3:])
        stat_df['Sim_Rho'] = sim_rho
        stat_df['Sim_Type'] = 'Selection'
        all_stat_dfs.append(stat_df)


#Now get the neutral data
#loop through each of the dataframes for the separate simulations
for curr_data in os.listdir(neutralDataDir):
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
    d_ratio_file = neutralDataDir + curr_data + "/linkage/d_ratio"
    stat_df = pd.read_pickle(d_ratio_file)
    stat_df['rep'] = int(rep[3:])
    stat_df['Sim_Rho'] = sim_rho
    stat_df['Sim_Type'] = 'Neutral'
    all_stat_dfs.append(stat_df)

all_stat_dfs = pd.concat(all_stat_dfs, ignore_index=True)
all_stat_dfs = all_stat_dfs[all_stat_dfs['Dist_X_Time'] <= DIST_TIME_MAX]

all_rho_bins = []

#Make the curves for the simulation types separately
for curr_type in all_stat_dfs['Sim_Type'].unique():
    type_stat_df = all_stat_dfs[all_stat_dfs['Sim_Type'] == curr_type]

    #Bin the points for each rho value
    for curr_rho in all_stat_dfs['Sim_Rho'].unique():
        curr_stat_df = type_stat_df[type_stat_df['Sim_Rho'] == curr_rho]
        if len(curr_stat_df) == 0:
            continue

        #Bin the d' ratios so they are easier to view on the plots
        binned_rat, binedges, bin_nums = binned_statistic(
            curr_stat_df['Dist_X_Time'].to_numpy().astype(float), 
            curr_stat_df['d_ratio'].to_numpy().astype(float), bins = 100)
        
        curr_rho_bins = pd.DataFrame({'Dist_X_Time': binedges[1:], 'D_Ratio': binned_rat})
        curr_rho_bins['Sim_Rho'] = curr_rho
        curr_rho_bins['Sim_Type'] = curr_type

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

############################# Plot the Data ###################################

print(all_rho_bins)
fig = sns.FacetGrid(data=all_rho_bins, col='Sim_Type', hue='Sim_float_rho', hue_kws = {'palette' : sns.color_palette("coolwarm", n_colors=len(all_rho_bins['Sim_float_rho'].unique()))})
fig.map(sns.lineplot, 'Dist_X_Time', 'D_Ratio',
             linewidth = 3,
             )

plt.tight_layout()
plt.legend()
plt.xlabel(r'Distance X Time (bp/generation)')
plt.ylabel("D\' Ratio")
print('Saving Figure')
plt.savefig(outDir + 'sim_fig_1C.png', dpi = 300)
plt.close()

############################# Plot the D' Distributions ###################################

example_rho = 0.0001
print(all_stat_dfs['Sim_Rho'].unique())
#print([type(x) for x in all_stat_dfs['Sim_Rho']])
my_example_df = all_stat_dfs[(all_stat_dfs['Sim_Rho'] == '1e-04')]
print(my_example_df.head())

fig = sns.FacetGrid(data=my_example_df, col='Sim_Type')
fig.map(sns.histplot, 'd_i', stat = 'density',
             linewidth = 3,
             )
plt.tight_layout()

plt.xlabel('D_i')
plt.ylabel("D\' Ratio")
print('Saving Figure')
plt.savefig(outDir + 'di_hist_faceted', dpi = 300)
plt.close()