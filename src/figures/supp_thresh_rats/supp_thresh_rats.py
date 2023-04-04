import sys
from tkinter import E
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
#for running on desktop
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plot_neher as plne
import autocorrelation as autocorr
import zaniniUtil as zu
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
outDir = "/net/feder/vol1/home/evromero/2021_hiv-rec/results/paper/supp_thresh_rats/"

# #For running on desktop
# dataDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/"
# vlDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/"
# outDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/paper/supp_thresh_rats/"

#Make the dataframe containg D' ratios
stat_df = zu.combine_drats(dataDir)
stat_df = zu.label_vl_drats(stat_df, vlDir)

#Filter the d' ratio dataframe
stat_df = stat_df[stat_df['Fragment'] != 'F5']
stat_df = stat_df[stat_df['Time_Diff'].gt(50)]
stat_df = stat_df[stat_df['Monotonic'] == True]
stat_df = stat_df[stat_df['d_i'] > THRESHOLD]

all_ests = []
all_group_fits = []
x_vals = stat_df['Dist_X_Time'].unique()

#Loop through the viral load thresholds
for curr_thresh in GROUP_THRESHOLD_LIST:
    print(curr_thresh)
    stat_df['High_VL'] = stat_df['Ave_VL'].gt(curr_thresh)

    #Estimate for the specific group
    for name, group in stat_df.groupby('High_VL'):
        #Name the group by whether it is high or low viral load
        if name:
            curr_name = 'High_VL'
        else: curr_name = 'Low_VL'

        #Get the current estimate
        lower_fit, upper_fit, estimate_df = plne.bootstrap_rho(group, NUM_BOOTSTRAPS)
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
        group_fits['Group'] = curr_name
        group_fits['Threshold'] = curr_thresh

        all_group_fits.append(group_fits)
        all_ests.append(estimate_df)

all_ests = pd.concat(all_ests, ignore_index=True)
all_group_fits = pd.concat(all_group_fits, ignore_index=True)

all_conf_ints = []
#Calculate the confidence intervals for the estimates with low + high VL
for name, group in all_ests.groupby(['Group', 'Threshold']):
    lower_conf = np.quantile(group['Estimated_Rho'], 0.025)
    upper_conf = np.quantile(group['Estimated_Rho'], 0.975)

    #append the confidence interval
    all_conf_ints.append([lower_conf, upper_conf, name[0], name[1]])
all_conf_ints = pd.DataFrame(all_conf_ints, 
    columns=['lower_conf', 'upper_conf', 'Group', 'Threshold'])

print(all_conf_ints)
all_conf_ints.to_csv(outDir + 'all_conf_ints' + str(NUM_BOOTSTRAPS)+'.csv')



########################### Plotting our results ##############################
params = {'font.size': 8, 'axes.labelsize': 8,'axes.titlesize':8,  'legend.fontsize': 8, 'xtick.labelsize': 8, 'ytick.labelsize': 8}
plt.rcParams.update(params)
linewidth = 1

fig, axs = plt.subplots(2,3, figsize = (6.5, 4), sharex = True, sharey = True)
myaxs = axs.flatten()

for i in range(len(GROUP_THRESHOLD_LIST)):
    curr_thresh = GROUP_THRESHOLD_LIST[i]
    curr_ax = myaxs[i]
    plot_df = all_group_fits[all_group_fits['Threshold'] == curr_thresh]
    sns.scatterplot(x ='bin_edges', y ='ratio_bins', hue = 'Group', palette = 'tab10', s = 5, ax = curr_ax, data = plot_df)

    #plot the confidence intervals separately so we don't mess up the legend
    high_df = plot_df[plot_df['Group'] == 'High_VL']
    low_df = plot_df[plot_df['Group'] == 'Low_VL']
    sns.lineplot(x ='bin_edges', y ='lower_conf', linewidth = linewidth, linestyle = 'dashed', color = 'saddlebrown', ax = curr_ax, data = high_df)
    sns.lineplot(x ='bin_edges', y ='mid_conf', linewidth = linewidth, color = 'saddlebrown', ax = curr_ax, data = high_df)
    sns.lineplot(x ='bin_edges', y ='high_conf', linewidth = linewidth, linestyle = 'dashed', color = 'saddlebrown', ax = curr_ax, data = high_df)
    sns.lineplot(x ='bin_edges', y ='lower_conf', linewidth = linewidth, linestyle = 'dashed', color = 'navy', ax = curr_ax, data = low_df)
    sns.lineplot(x ='bin_edges', y ='mid_conf', linewidth = linewidth, color = 'navy', ax = curr_ax, data = low_df)
    sns.lineplot(x ='bin_edges', y ='high_conf', linewidth = linewidth, linestyle = 'dashed', color = 'navy', ax = curr_ax, data = low_df)
    curr_ax.set_xlabel(r'Distance $\cdot$' + "Time" + r'(bp $\cdot$ generation)')
    curr_ax.set_ylabel("D' Ratio")
    curr_ax.set_title('Threshold = ' + str(curr_thresh))

    if i != 0:
        curr_ax.get_legend().remove()
    else:
        new_labels = ['Low VL', 'High VL']
        for t, l in zip(curr_ax.get_legend().texts, new_labels):
            t.set_text(l)
        for curr_dot in curr_ax.legend_.legendHandles:
            curr_dot._sizes = [5]

# darker = ['navy', 'saddlebrown']
# #sns.set(rc = {'font.size': 8,'figure.figsize':(6.5, 4) })
# linewidth = 1
# myplot = sns.FacetGrid(all_group_fits, col = 'Threshold', col_wrap = 3, height = 2, aspect = 1, sharey = True, sharex = True)
# myplot.map_dataframe(sns.scatterplot, x ='bin_edges', y ='ratio_bins', hue = 'Group', palette = 'tab10', s = 3)
# plt.legend()
# myplot.map_dataframe(sns.lineplot, x ='bin_edges', y ='lower_conf', linewidth = linewidth, linestyle = 'dashed', hue = 'Group', palette = sns.color_palette(darker))
# myplot.map_dataframe(sns.lineplot, x ='bin_edges', y ='mid_conf', linewidth = linewidth, hue = 'Group', palette = sns.color_palette(darker))
# myplot.map_dataframe(sns.lineplot, x ='bin_edges', y ='high_conf', linewidth = linewidth, linestyle = 'dashed', hue = 'Group', palette = sns.color_palette(darker))
# myplot.set_titles('{col_name} copies/mL')
# myplot.set_xlabels(r'Distance $\cdot$' + "Time \n" + r'(bp $\cdot$ generation)')
# #for some reason the y label latex only works if you make the string type ""
# myplot.set_ylabels("D' Ratio")
plt.tight_layout()
plt.savefig(outDir + "fits_sep_groups" + str(NUM_BOOTSTRAPS)+".jpg", dpi = 300)
plt.close()