import sys
# sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
#for running on desktop
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic
import plot_neher as plne
import autocorrelation as autocorr
from scipy import optimize


#This script makes figure 1 for my general exam proposal

#This script just makes a basic plot of autoccorelation estimates for simulated data against
#the actual values in the simulation
#from this analysis, it looks like fitting on the sampled data doesn't quite work
THRESHOLD = 0.2
SEG_LOCI = 100

#For running on cluster
dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_02_24/'
outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_02_24/'

dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_02_24/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_02_24/05-06-2022/'

estimate_df = [] 

for curr_data in os.listdir(dataDir):
    #only get the data directories, not hidden files
    if curr_data[0] == '.':
        continue
    run_info = curr_data.split('_')
    sim_rho = run_info[1]
    sim_rho = sim_rho[3:]
    rep = run_info[-1]

    #make a place to store our output
    currOut = outDir + curr_data
    if not os.path.exists(currOut):
        os.mkdir(currOut)

    linkage_file = dataDir + curr_data + "/linkage/r2_and_D"
    d_ratio_file = dataDir + curr_data + "/linkage/d_ratio"

    if os.path.exists(d_ratio_file):
        stat_df = pd.read_pickle(d_ratio_file)
    else:
        stat_df = autocorr.calculate_d_ratios(linkage_file, THRESHOLD)
        stat_df.to_pickle(d_ratio_file)
    
    #make a list of all the segregating loci
    all_seg_loc_1 = set(stat_df["Locus_1"].unique())
    all_seg_loc_2 = set(stat_df["Locus_2"].unique())
    all_seg_loc = all_seg_loc_1.union(all_seg_loc_2)
    all_seg_loc = np.array(list(all_seg_loc))

    #sample a given number of segregating loci
    seg_loc_sample = np.random.choice(all_seg_loc, size = SEG_LOCI, replace = False)
    seg_loc_sample = set(seg_loc_sample)

    #get only the autocorrelation of the chosen loci
    stat_df_sample = stat_df[stat_df["Locus_1"].isin(seg_loc_sample)]
    stat_df_sample = stat_df_sample[stat_df_sample["Locus_2"].isin(seg_loc_sample)]

    x_vals = stat_df_sample['Dist_X_Time'].unique()
    coeffs, fit_dat = optimize.curve_fit(plne.neher_leitner, stat_df_sample['Dist_X_Time'], stat_df_sample['d_ratio'], p0 = [0, 0.26, .0000439])
    fit_vals = [plne.neher_leitner(x, coeffs[0], coeffs[1], coeffs[2]) for x in x_vals]
    # sns.scatterplot(x = 'Dist_X_Time', y = 'd_ratio', data = stat_df_sample)
    # sns.lineplot(x = 'Dist_X_Time', y = 'd_ratio', data = stat_df_sample, estimator = np.mean)
    # sns.lineplot(x = 'Dist_X_Time', y = 'd_ratio', data = stat_df_sample, estimator = np.median)
    # plt.savefig(currOut + "/auto_plot.jpg")
    # plt.close()

    estimate_df.append([coeffs[0], coeffs[1], coeffs[2], coeffs[1] * coeffs[2], curr_data, sim_rho, 'curve', curr_data])  

estimate_df = pd.DataFrame(estimate_df, columns=["C0", "C1", "C2", "Est_Rho", 'Dataset', 'Sim_Rho', 'style', 'data'] )

############################# Plotting Estimate Accuracy ######################
# #Plot our estimates against each other 
#make the rho values ints rather than strings
rhoDict = {"0.001" : 0.001,
            "1e-04" : 0.0001,
            "2e-04" : 0.0002,
            "1e-05" : 0.00001,
            "2e-05" : 0.00002,
            "2e-06" : 0.000002}

intRhoList = []
for entry in estimate_df['Sim_Rho']:
    intRhoList.append(rhoDict[entry])
estimate_df['Sim_int_rho'] = intRhoList

fig, ax = plt.subplots(2)
sns.stripplot(x = 'Sim_Rho', y = 'Est_Rho', data = estimate_df, jitter = True, color = 'k', s = 4, alpha = 0.5,
    order = ["0.001", "2e-04", "1e-04", "2e-05", "1e-05", "2e-06"], ax = ax[0])

# distance across the "X" or "Y" stipplot column to span, in this case 40%
label_width = 0.4

for tick, text in zip(ax[0].get_xticks(), ax[0].get_xticklabels()):
    sample_name = text.get_text()  # "X" or "Y"

    # calculate the median value for all replicates of either X or Y
    rho_val = rhoDict[sample_name]

    # plot horizontal lines across the column, centered on the tick
    ax[0].plot([tick-label_width/2, tick+label_width/2], [rho_val, rho_val],
            lw=2, color='k')

ax[0].legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
ax[0].set_xlabel(r'Simulation Value of $\rho$')
ax[0].set_ylabel(r'Estimated Value of $\rho$')
ax[0].set_ylim(0.000001, 0.01)
ax[0].set_yscale('log')

#I'm going to plot d ratios over the length of the paired end reads

dataDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/"
vlDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/"
outDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/zanini/05-23-2022/"

#loop through the files and plot the numbers of segregating loci
all_d_rats = []

#loop through all of the directories with linkage files
for curr_dir in os.listdir(dataDir):
    if curr_dir[0] == '.':
        continue
    for curr_file in os.listdir(dataDir + curr_dir):
        #get the participant and timepoint
        curr_par = curr_dir.split('_')[0]
        curr_frag = curr_dir.split('_')[1]

        #get the timepoint dictionary for the current participant
        curr_seg = pd.read_pickle(dataDir + curr_dir + "/linkage/d_ratio")
        curr_seg['Participant'] = curr_par
        curr_seg['Fragment'] = curr_frag
        all_d_rats.append(curr_seg)

#put all the ratios together
all_d_rats = pd.concat(all_d_rats, ignore_index= True)
all_d_rats['d_i'] = all_d_rats['d_i'].to_numpy().astype(float)
all_d_rats['d_ratio'] = all_d_rats['d_ratio'].to_numpy().astype(float)
all_d_rats.sort_values('Dist_X_Time', inplace= True)
all_d_rats['d_i_1'] = all_d_rats['d_i'] * np.exp(np.negative(all_d_rats['d_ratio']))

all_d_rats = all_d_rats[all_d_rats['Participant'] != 'p10']
print(all_d_rats)
all_bins = []
for curr_par in all_d_rats['Participant'].unique():
    curr_d_rats = all_d_rats[all_d_rats['Participant'] == curr_par]
    #don't include samples from time points too close together (we previously noted noise issue)
    #with simulated data
    curr_d_rats = curr_d_rats[curr_d_rats['Time_Diff'].gt(50)]
    dist_time_vals = curr_d_rats['Dist_X_Time'].to_numpy().astype(float)
    d_i_vals = curr_d_rats['d_i']
    d_i_1_vals = curr_d_rats['d_i_1']
    all_d_vals = np.concatenate((d_i_1_vals, d_i_vals))
    all_dist_time = np.concatenate((dist_time_vals, dist_time_vals))
    dist_time_vals = curr_d_rats['Dist_X_Time'].to_numpy().astype(float)
    rolling_ave, bin_edges, binnumber = binned_statistic(all_dist_time, all_d_vals, bins = 25, statistic = 'mean')

    #plot things at the center of bins
    bin_edges = pd.Series(bin_edges)
    bin_edges = bin_edges.rolling(2).mean()
    bin_edges = bin_edges.to_numpy()[1:]
    par_bin_data = pd.DataFrame(zip(rolling_ave, bin_edges), columns = ['average', 'edges'])
    par_bin_data['Participant'] = curr_par
    all_bins.append(par_bin_data)
all_bins = pd.concat(all_bins, ignore_index=True)


# sns.scatterplot(x = 'Dist_X_Time', y = 'd_i', data = all_d_rats, alpha = 0.007, ax = ax[1], hue = 'Participant')
sns.lineplot(y = 'average', x = 'edges', data = all_bins, linewidth = 4, ax = ax[1], hue = 'Participant')
ax[1].set_xlabel('Distance x Time (bp x generations)' )
ax[1].set_ylabel('D\' Value')
plt.tight_layout()
plt.legend(bbox_to_anchor=(1.3,0.3), loc="lower right")
# plt.subplots_adjust(left=0.1, bottom=0.1, right=0.75)

plt.savefig(outDir + "figure_1_mean.jpg", dpi = 300)
plt.close()