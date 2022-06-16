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
import zaniniUtil as zu
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

#make the dataframe containg D' ratios
stat_df = zu.combine_drats(dataDir)
stat_df = zu.label_vl_drats(stat_df, vlDir)

#filter the d' ratio dataframe
stat_df = stat_df[stat_df['Fragment'] != 'F5']
stat_df = stat_df[stat_df['Time_Diff'].gt(50)]
stat_df = stat_df[stat_df['Participant'] != 'p7']
stat_df = stat_df[stat_df['Participant'] != 'p4']
stat_df = stat_df[stat_df['Participant'] != 'p10']

lower_fit, upper_fit, estimate_df = plne.bootstrap_rho(stat_df, 1000)
print(lower_fit)

final_estimates = []

############# Plotting our results###############################
x_vals = stat_df['Dist_X_Time'].unique()
binned_rat, binedges, bin_nums = stats.binned_statistic(stat_df['Dist_X_Time'].to_numpy(), stat_df['d_ratio'].to_numpy(), bins = 100)


#Fit the data and add our estimate to the dataframe
coeffs, fit_dat = optimize.curve_fit(plne.neher_leitner, stat_df['Dist_X_Time'], stat_df['d_ratio'], p0 = [0, 0.26, .0000439])
lower_fit, upper_fit, all_boots2 = plne.bootstrap_rho(stat_df)
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