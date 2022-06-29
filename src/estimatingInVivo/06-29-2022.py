import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
#for running on desktop
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plot_neher as plne
import zaniniUtil as zu
from scipy import stats

# In this script I am repeating the viral load analysis, but I am resampling
# from the patients such that I leave one out in each sample
THRESHOLD = 0.2
DIST_TIME_MAX = 50000
GROUP_THRESHOLD_LIST = [300, 400, 500, 600, 700, 800, 900, 1000]
NUM_BOOTSTRAPS = 1000

#For running on desktop
dataDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/"
cd4Dir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini/cd4Counts/"
vlDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/"
outDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/zanini/06-29-2022/"

#Make the dataframe containg D' ratios
stat_df = zu.combine_drats(dataDir)
stat_df = zu.label_vl_drats(stat_df, cd4Dir, cd4_label = True)

#filter the d' ratio dataframe
stat_df = stat_df[stat_df['Fragment'] != 'F5']
stat_df = stat_df[stat_df['Time_Diff'].gt(50)]
stat_df = stat_df[stat_df['Participant'] != 'p10']
stat_df = stat_df[stat_df['Monotonic'] == True]

stat_df = zu.label_vl_drats(stat_df, vlDir, cd4_label = False)
stat_df = stat_df[stat_df['Monotonic'] == True]

all_par_ests = []
all_group_fits = []
x_vals = stat_df['Dist_X_Time'].unique()

################### Plot the CD4+ Counts ######################################
sns.histplot(x = 'Ave_CD4', data = stat_df)
plt.savefig(outDir + "CD4_hist.jpg")
plt.close()

sns.scatterplot(x = 'Ave_CD4', y = 'Ave_VL', data = stat_df)
plt.yscale('log')
plt.savefig(outDir + "CD4_VL_scatter.jpg")


#Estimate rates specifically excluding each individual
for curr_thresh in GROUP_THRESHOLD_LIST:
    #Get the dataframe for everyone except the current participant
    stat_df['High_CD4'] = stat_df['Ave_CD4'].gt(curr_thresh)

    #Estimate for the specific group
    for name, group in stat_df.groupby('High_CD4'):
        #Name the group by whether it is high or low viral load
        if name:
            curr_name = 'High_CD4'
        else: curr_name = 'Low_CD4'


        #Get the current estimate
        lower_fit, upper_fit, estimate_df = plne.bootstrap_rho(group, NUM_BOOTSTRAPS)
        estimate_df['Threshold'] = curr_thresh
        estimate_df['Group'] = curr_name

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
        group_fits['Threshold'] = curr_thresh
        group_fits['Group'] = curr_name

        all_group_fits.append(group_fits)
        all_par_ests.append(estimate_df)

all_par_ests = pd.concat(all_par_ests, ignore_index=True)
all_group_fits = pd.concat(all_group_fits, ignore_index=True)

all_conf_ints = []
#Calculate the confidence intervals for the estimates with low + high CD4
for name, group in all_par_ests.groupby(['Group', 'Threshold']):
    lower_conf = np.quantile(group['Estimated_Rho'], 0.025)
    upper_conf = np.quantile(group['Estimated_Rho'], 0.975)

    #append the confidence interval
    all_conf_ints.append([lower_conf, upper_conf, name[0], name[1]])
all_conf_ints = pd.DataFrame(all_conf_ints, 
    columns=['lower_conf', 'upper_conf', 'Group', 'Threshold'])

print(all_conf_ints)




########################### Plotting our results ##############################
myplot = sns.FacetGrid(all_group_fits, col = 'Group', row = 'Threshold')
myplot.map_dataframe(sns.scatterplot, x ='bin_edges', y ='ratio_bins')
myplot.map_dataframe(sns.lineplot, x ='bin_edges', y ='lower_conf', color = 'gray', linestyle = 'dashed')
myplot.map_dataframe(sns.lineplot, x ='bin_edges', y ='mid_conf', color = 'k')
myplot.map_dataframe(sns.lineplot, x ='bin_edges', y ='high_conf', color = 'gray', linestyle = 'dashed')
plt.xlabel(r'Distance X Time (bp/generation)')
plt.ylabel(r'D\' Ratio')
plt.savefig(outDir + "thresh_fits_sep_groups.jpg")
plt.close()

################### Estimate Facet by Excl_Par ################################
myplot = sns.FacetGrid(all_par_ests, col = 'Threshold', col_wrap = 2)
myplot.map_dataframe(sns.boxplot, x ='Group', y ='Estimated_Rho')
plt.xlabel(r'Group')
plt.ylabel(r'Estimated Value of $\rho$')
plt.tight_layout()
plt.savefig(outDir + "thresh_estimates.jpg")
plt.close()

