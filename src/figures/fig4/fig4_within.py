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

# In this script I am estimating the recombination rate for p1 at low and high
# viral load timepoints.
THRESHOLD = 0.2
DIST_TIME_MAX = 50000
GROUP_THRESH = 50000
NUM_BOOTSTRAPS = 1000


#For running on Cluster
dataDir = "/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/"
vlDir = "/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/"
outDir = "/net/feder/vol1/home/evromero/2021_hiv-rec/results/paper/fig4/"

#For running on desktop
dataDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/"
vlDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/"
outDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/paper/fig4/"

#Make the dataframe containg D' ratios
all_stat_df = zu.combine_drats(dataDir)
all_stat_df = zu.label_vl_drats(all_stat_df, vlDir)
all_stat_df = all_stat_df[all_stat_df['Monotonic'] == True]
print(np.unique(all_stat_df['Participant'], return_counts = True))


#Filter the d' ratio dataframe
all_stat_df = all_stat_df[all_stat_df['Time_Diff'].gt(50)]
par_list = ['p3']

for curr_par in par_list:
    stat_df = all_stat_df[all_stat_df['Participant'] == curr_par]  
    stat_df['High_VL'] = stat_df['Ave_VL'].gt(GROUP_THRESH)

    print(stat_df['Participant'].unique())
    print(stat_df.columns)

    all_par_ests = []
    all_group_fits = []
    x_vals = stat_df['Dist_X_Time'].unique()


    #Estimate Rates Specifically for each individual
    for name, group in stat_df.groupby('High_VL'):
        #Name the group by whether it is high or low viral load
        if name:
            curr_name = 'High_VL'
        else: curr_name = 'Low_VL'
        
        #Get estimates for the group
        lower_fit, upper_fit, estimate_df = plne.bootstrap_rho(group, NUM_BOOTSTRAPS)
        if name:
            estimate_df['Group'] = curr_name
        else: estimate_df['Group'] = curr_name

        #Bin the d' ratios so they are easier to view on the plots
        binned_rat, binedges, bin_nums = stats.binned_statistic(
            group['Dist_X_Time'].to_numpy(), 
            group['d_ratio'].to_numpy(), bins = 100)
        
        #troubleshoot with plot
        sns.scatterplot(x = binedges[1:], y = binned_rat, label = curr_name)
        plt.xlabel(r'Distance X Time (bp/generation)')
        plt.ylabel(r'D\' Ratio')

        #Get the mean value of the estimates
        mid_fit = np.quantile(estimate_df['Estimated_Rho'], 0.5)
        mid_fit = estimate_df.iloc[(estimate_df['Estimated_Rho']-mid_fit).abs().argsort()[:2]]

        #Get the low and high values for the confidence interval
        fit_vals_low = [plne.neher_leitner(x, lower_fit['c0'].to_numpy()[0], lower_fit['c1'].to_numpy()[0], lower_fit['c2'].to_numpy()[0]) for x in binedges]
        fit_vals_high = [plne.neher_leitner(x, upper_fit['c0'].to_numpy()[0], upper_fit['c1'].to_numpy()[0], upper_fit['c2'].to_numpy()[0]) for x in binedges]
        fit_vals_mid = [plne.neher_leitner(x, mid_fit['c0'].to_numpy()[0], mid_fit['c1'].to_numpy()[0], mid_fit['c2'].to_numpy()[0]) for x in binedges]

        #Gather all of the fits for the participant
        group_fits = pd.DataFrame(zip(binned_rat, binedges, fit_vals_low, fit_vals_mid, fit_vals_high),
                        columns=['ratio_bins', 'bin_edges', 'lower_conf', 'mid_conf', 'high_conf'])
        group_fits['Group'] = curr_name

        all_group_fits.append(group_fits)
        all_par_ests.append(estimate_df)

        
    all_par_ests = pd.concat(all_par_ests, ignore_index=True)
    all_group_fits = pd.concat(all_group_fits, ignore_index=True)


    ########################### Plotting our results ##############################
    myplot = sns.FacetGrid(all_group_fits, col = 'Group')
    myplot.map_dataframe(sns.scatterplot, x ='bin_edges', y ='ratio_bins')
    myplot.map_dataframe(sns.lineplot, x ='bin_edges', y ='lower_conf', color = 'gray', linestyle = 'dashed')
    myplot.map_dataframe(sns.lineplot, x ='bin_edges', y ='mid_conf', color = 'k')
    myplot.map_dataframe(sns.lineplot, x ='bin_edges', y ='high_conf', color = 'gray', linestyle = 'dashed')
    plt.xlabel(r'Distance X Time (bp/generation)')
    plt.ylabel(r'D\' Ratio')
    plt.savefig(outDir + curr_par + "_fits_by_vl.jpg")
    plt.close()


    ################################ Estimate Plot ################################
    sns.set(rc={'figure.figsize':(15,7.5)}, font_scale = 2)
    sns.set_palette("tab10")
    ax = sns.boxplot(x = 'Group', y = 'Estimated_Rho', data = all_par_ests)
    ax.axhline(0.000008, linestyle = 'dashed', color = 'tab:green')
    ax.axhline(0.000014, color = 'tab:green')
    ax.axhline(0.00002, linestyle = 'dashed', color = 'tab:green')
    ax.set_ylim(0, 0.0002)
    # distance across the "X" or "Y" stipplot column to span, in this case 40%
    label_width = 0.4

    plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
    plt.xlabel(r'Group')
    plt.ylabel(r'Estimated Value of $\rho$')
    plt.tight_layout()
    plt.savefig(outDir + curr_par + "_estimates_by_group.jpg")
    plt.close()
