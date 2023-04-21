import stat
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
from matplotlib import rcParams

# In this script I am repeating the viral load analysis, but I am resampling
# from the patients such that I leave one out in each sample
THRESHOLD = 0.2
DIST_TIME_MAX = 50000
GROUP_THRESHOLD_LIST = [7500, 10000, 25000, 50000, 100000, 200000]
NUM_BOOTSTRAPS = 1000

#For running on Cluster
dataDir = "/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/"
vlDir = "/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/"
outDir = "/net/feder/vol1/home/evromero/2021_hiv-rec/results/paper/supp_loo/"

#For running on desktop
dataDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/"
vlDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/"
outDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/paper/supp_loo/"

#Make the dataframe containg D' ratios
stat_df = zu.combine_drats(dataDir)
stat_df = zu.label_vl_drats(stat_df, vlDir)

#Filter the d' ratio dataframe
stat_df = stat_df[stat_df['Fragment'] != 'F5']
# stat_df = stat_df[stat_df['Participant'] != 'p4']
stat_df = stat_df[stat_df['Time_Diff'].gt(50)]
stat_df = stat_df[stat_df['Monotonic'] == True]
stat_df = stat_df[stat_df['d_i'] > THRESHOLD]

all_par_ests = []
all_group_fits = []
x_vals = stat_df['Dist_X_Time'].unique()

#Loop through the viral load thresholds
for curr_thresh in GROUP_THRESHOLD_LIST:
    print(curr_thresh)
    stat_df['High_VL'] = stat_df['Ave_VL'].gt(curr_thresh)

    #Estimate rates specifically excluding each individual
    for curr_par in stat_df['Participant'].unique():
        #Get the dataframe for everyone except the current participant
        curr_par_df = stat_df[stat_df['Participant'] != curr_par]

        #Estimate for the specific group
        for name, group in curr_par_df.groupby('High_VL'):
            #Name the group by whether it is high or low viral load
            if name:
                curr_name = 'High_VL'
            else: curr_name = 'Low_VL'

            #Get the current estimate
            lower_fit, upper_fit, estimate_df = plne.bootstrap_rho(group, NUM_BOOTSTRAPS)
            estimate_df['excl_par'] = curr_par
            estimate_df['Group'] = curr_name
            estimate_df['Threshold'] = curr_thresh

            #Bin the d' ratios so they are easier to view on the plots
            binned_rat, binedges, bin_nums = stats.binned_statistic(
                group['Dist_X_Time'].to_numpy(), 
                group['d_ratio'].to_numpy(), bins = 100)

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
            group_fits['excl_par'] = curr_par
            group_fits['Group'] = curr_name
            group_fits['Threshold'] = curr_thresh

            all_group_fits.append(group_fits)
            all_par_ests.append(estimate_df)

all_par_ests = pd.concat(all_par_ests, ignore_index=True)
all_group_fits = pd.concat(all_group_fits, ignore_index=True)

all_conf_ints = []
#Calculate the confidence intervals for the estimates with low + high VL
for name, group in all_par_ests.groupby(['Group', 'excl_par', 'Threshold']):
    lower_conf = np.quantile(group['Estimated_Rho'], 0.025)
    upper_conf = np.quantile(group['Estimated_Rho'], 0.975)

    #append the confidence interval
    all_conf_ints.append([lower_conf, upper_conf, name[0], name[1], name[2]])
all_conf_ints = pd.DataFrame(all_conf_ints, 
    columns=['lower_conf', 'upper_conf', 'Group', 'excl_par', 'Threshold'])

print(all_conf_ints)
all_conf_ints.to_csv(outDir + 'all_conf_ints' + str(NUM_BOOTSTRAPS)+'.csv')



########################### Plotting our results ##############################
linewidth = 1
myplot = sns.FacetGrid(all_group_fits, col = 'Threshold', row = 'excl_par')
myplot.map_dataframe(sns.scatterplot, x ='bin_edges', y ='ratio_bins', hue = 'Group', palette = 'tab10')
myplot.map_dataframe(sns.lineplot, x ='bin_edges', y ='lower_conf',linestyle = 'dashed', hue = 'Group', palette = 'tab10')
myplot.map_dataframe(sns.lineplot, x ='bin_edges', y ='mid_conf', hue = 'Group', palette = 'tab10')
myplot.map_dataframe(sns.lineplot, x ='bin_edges', y ='high_conf', linestyle = 'dashed', hue = 'Group', palette = 'tab10')
plt.xlabel(r'Distance X Time (bp/generation)')
plt.ylabel(r'D\' Ratio')
plt.savefig(outDir + "loo_fits_sep_groups.jpg")
plt.close()

################### Estimate Facet by Excl_Par ################################
params = {'axes.labelsize': 8,'axes.titlesize':8,  'legend.fontsize': 8, 'xtick.labelsize': 8, 'ytick.labelsize': 8}
plt.rcParams.update(params)

#We need to use facetgrid params to adjust the size of the plots
#Because outherwise seaborn just writes over the figure size
myplot = sns.FacetGrid(all_par_ests, col = 'excl_par', col_order = ['p1', 'p2', 'p3','p4', 'p5', 'p6', 'p7', 'p8', 'p9', 'p10', 'p11'], col_wrap = 3,
    height = 2, aspect = 1)
myplot.map_dataframe(sns.boxplot, x ='Threshold', y ='Estimated_Rho', hue = 'Group', palette = 'tab10', fliersize = 0, linewidth = linewidth)


# axd['right'].legend(bbox_to_anchor=(1, 0.5))
new_labels = ['Low VL', 'High VL']
for t, l in zip(plt.legend().texts, new_labels):
    t.set_text(l)

#make the layout tight before we add any labels or legends that would push the
#panels further apart.
plt.tight_layout()

myplot.set_ylabels("Estimated Recombination \n" + r'Rate ($\hat{\rho}$)')
myplot.set_xlabels(r'Group Viral Load Threshold' + "\n (copies/mL)")
myplot.set_titles('Excluding {col_name}')


labeledaxes = [8,9,10]
for curr_ind in labeledaxes:
    labels = myplot.axes[curr_ind].get_xticklabels()
    myplot.axes[curr_ind].set_xticklabels(labels, rotation=45) # set new labels

myplot.axes[10].legend(loc='center left', bbox_to_anchor=(1,0.5))
new_labels = ['Low VL', 'High VL']
for t, l in zip(myplot.axes[10].get_legend().texts, new_labels):
    t.set_text(l)

plt.ticklabel_format(style = 'sci', axis = 'y', scilimits = (0,0))
plt.ylim(0, 0.0002)
plt.savefig(outDir + "loo_estimates_by_excl_" + str(NUM_BOOTSTRAPS) + ".jpg", dpi = 300)
plt.close()