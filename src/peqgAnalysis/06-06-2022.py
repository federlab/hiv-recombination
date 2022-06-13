import sys
from tkinter import E
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
#for running on desktop
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import os
import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plot_neher as plne
import autocorrelation as autocorr
from scipy import optimize
from scipy import stats

#In this script I am going to plot the average D' ratio as a function of distance and time

#I tried measuring the autocorrelation of the D statistic, but it looks like we
#are getting a lot of noise. So I am going to try setting up an an initial 
#thresholding value to only run tests after high linkage is initially seen.
THRESHOLD = 0.2
DIST_TIME_MAX = 50000
GROUP_THRESH = 50000


#For running on Cluster
dataDir = "/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/"
vlDir = "/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/"
outDir = "/net/feder/vol1/home/evromero/2021_hiv-rec/results/zanini/poster_peqg/"

#For running on desktop
dataDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/"
vlDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/"
outDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/zanini/poster_peqg/"

stat_df = []

for curr_data in os.listdir(dataDir):
    if curr_data[0] == '.':
        continue
    # #only get fragment 5 for now
    # if curr_data.split('_')[1] != 'F5':
    #     continue

    #
    run_info = curr_data.split('_')
    curr_par = run_info[0]
    curr_frag = run_info[1]

    linkage_file = dataDir + curr_data + "/linkage/r2_and_D"
    d_ratio_file = dataDir + curr_data + "/linkage/d_ratio"

    curr_stat_df = pd.read_pickle(d_ratio_file)
    curr_stat_df['Participant'] = curr_par
    curr_stat_df['Fragment'] = curr_frag
    stat_df.append(curr_stat_df)


#put all the ratios together
stat_df = pd.concat(stat_df, ignore_index= True)
print(stat_df.columns)
stat_df['d_i'] = stat_df['d_i'].to_numpy().astype(float)
stat_df['d_ratio'] = stat_df['d_ratio'].to_numpy().astype(float)
stat_df['Dist_X_Time'] = stat_df['Dist_X_Time'].to_numpy().astype(float)
stat_df['d_i_1'] = stat_df['d_i'] * np.exp(np.negative(stat_df['d_ratio']))
stat_df = stat_df[stat_df['Fragment'] != 'F5']
stat_df = stat_df[stat_df['Time_Diff'].gt(50)]
stat_df = stat_df[stat_df['Participant'] != 'p7']
stat_df = stat_df[stat_df['Participant'] != 'p4']
stat_df = stat_df[stat_df['Participant'] != 'p10']
# stat_df = stat_df[stat_df['Fragment'] == 'F5']
stat_df = stat_df[stat_df['Dist_X_Time'].between(0, DIST_TIME_MAX)]
stat_df['Dist'] = stat_df['Locus_2'] - stat_df['Locus_1']
print("Max dist is: " + str(max(stat_df['Dist'])))

labeled_rats = []
#label each participants ratios with the corresponding viral loads
for curr_file in os.listdir(vlDir):
    dict_from_csv = {}
    if curr_file[0] == '.':
        continue

    #make a dictionary of timepoints and their viral loads
    with open(vlDir + curr_file, mode='r') as inp:
        reader = csv.reader(inp, delimiter= '\t')
        next(reader)
        dict_from_csv = {float(rows[0]):float(rows[1]) for rows in reader}
    
    #get the ratios for the participant
    participant = curr_file.split('.')[0]
    participant = participant.split('_')[1]
    curr_d_rats = stat_df[stat_df['Participant'] == participant]
    curr_d_rats['Day_1'] = curr_d_rats['Time_1'] * 2
    curr_d_rats['Day_2'] = curr_d_rats['Time_2'] * 2

    #label the ratios
    curr_d_rats = curr_d_rats[curr_d_rats['Day_1'].isin(dict_from_csv.keys())]
    curr_d_rats = curr_d_rats[curr_d_rats['Day_2'].isin(dict_from_csv.keys())]

    curr_d_rats['VL_1'] = curr_d_rats['Day_1'].map(lambda x: dict_from_csv[x])
    curr_d_rats['VL_2'] = curr_d_rats['Day_2'].map(lambda x: dict_from_csv[x])
    labeled_rats.append(curr_d_rats)

labeled_rats = pd.concat(labeled_rats, ignore_index= True)
labeled_rats['Ave_VL'] = labeled_rats[['VL_1', 'VL_2']].mean(axis=1)
stat_df = labeled_rats


def bootstrap_rho(d_ratio_df):
    """Takes a d_ratio dataframe and resamples from it to bootstrap a rho estimate"""
    #make a list of all the segregating loci
    all_seg_loc_1 = set(d_ratio_df["Locus_1"].unique())
    all_seg_loc_2 = set(d_ratio_df["Locus_2"].unique())
    all_seg_loc = all_seg_loc_1.union(all_seg_loc_2)
    all_seg_loc = np.array(list(all_seg_loc))
    num_seg = len(all_seg_loc)

    estimates = []

    #sample with replacement 10k times
    for i in range(1000):
        #sample a given number of segregating loci
        seg_loc_sample = np.random.choice(all_seg_loc, size = num_seg, replace = True)
        seg_loc_sample = set(seg_loc_sample)

        #get only the autocorrelation of the chosen loci
        stat_df_sample = d_ratio_df[d_ratio_df["Locus_1"].isin(seg_loc_sample)]
        stat_df_sample = stat_df_sample[stat_df_sample["Locus_2"].isin(seg_loc_sample)]

        coeffs, fit_dat = optimize.curve_fit(plne.neher_leitner, stat_df_sample['Dist_X_Time'], stat_df_sample['d_ratio'], p0 = [0, 0.26, .0000439])
        estimates.append([coeffs[0], coeffs[1], coeffs[2], coeffs[1] * coeffs[2]])
    estimate_df = pd.DataFrame(estimates, columns = ['c0', 'c1', 'c2', 'Estimated_Rho'])
    conf_int = (np.quantile(estimate_df['Estimated_Rho'], 0.025), np.quantile(estimate_df['Estimated_Rho'], 0.975))
    print(np.quantile(estimate_df['Estimated_Rho'], 0.5))
    #get the row at the given quantile
    lower_fit = estimate_df.iloc[(estimate_df['Estimated_Rho']-conf_int[0]).abs().argsort()[:2]]
    lower_fit = lower_fit.head(1)
    upper_fit = estimate_df.iloc[(estimate_df['Estimated_Rho']-conf_int[1]).abs().argsort()[:2]]
    upper_fit = upper_fit.head(1)
    return lower_fit, upper_fit, estimate_df

final_estimates = []

############# Plotting our results###############################
x_vals = stat_df['Dist_X_Time'].unique()
binned_rat, binedges, bin_nums = stats.binned_statistic(stat_df['Dist_X_Time'].to_numpy(), stat_df['d_ratio'].to_numpy(), bins = 100)


#Fit the data and add our estimate to the dataframe
coeffs, fit_dat = optimize.curve_fit(plne.neher_leitner, stat_df['Dist_X_Time'], stat_df['d_ratio'], p0 = [0, 0.26, .0000439])
lower_fit, upper_fit, all_boots2 = bootstrap_rho(stat_df)
fit_vals_low = [plne.neher_leitner(x, lower_fit['c0'].to_numpy()[0], lower_fit['c1'].to_numpy()[0], lower_fit['c2'].to_numpy()[0]) for x in x_vals]
fit_vals_high = [plne.neher_leitner(x, upper_fit['c0'].to_numpy()[0], upper_fit['c1'].to_numpy()[0], upper_fit['c2'].to_numpy()[0]) for x in x_vals]
fit_vals = [plne.neher_leitner(x, coeffs[0], coeffs[1], coeffs[2]) for x in x_vals]
sns.lineplot(x = x_vals, y = fit_vals, color = 'tab:blue', linestyle = 'dashed', linewidth = 3)
sns.lineplot(x = binedges[:-1], y = binned_rat, color = 'tab:blue')
sns.lineplot(x = x_vals, y = fit_vals_low, color = 'black')
sns.lineplot(x = x_vals, y = fit_vals_high, color = 'black')
final_estimates.append([coeffs[1] * coeffs[2], lower_fit['Estimated_Rho'].to_numpy()[0], upper_fit['Estimated_Rho'].to_numpy()[0]])
final_estimates = pd.DataFrame(final_estimates, columns = ['Estimate', 'Lower', 'Upper'])
print(final_estimates)
plt.ylabel('D\' ratio')
plt.xlabel('Distance X Time (bp x generations)')
plt.savefig(outDir + "auto_plot_binned_vl_except_frag_5_" + str(DIST_TIME_MAX) +".jpg", dpi = 300)
plt.close()
