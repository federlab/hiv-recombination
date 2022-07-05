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
from matplotlib import rcParams

LOWER_THRESHOLD = 50000
UPPER_THRESHOLD = 100000
NUM_BOOTSTRAPS = 10

#In this script, I am splitting the data into three groups and estimating the 
#viral load for each group.

#For running on desktop
dataDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/"
vlDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/"
outDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/zanini/06-29-2022/"

#Make the dataframe containg D' ratios
stat_df = zu.combine_drats(dataDir)
stat_df = zu.label_vl_drats(stat_df, vlDir)

stat_df = stat_df[stat_df['Fragment'] != 'F5']
stat_df = stat_df[stat_df['Participant'] != 'p4']
stat_df = stat_df[stat_df['Time_Diff'].gt(50)]

#Only get values for which the viral load doesn't leave the boundaries of the 
#initial and final measurments
stat_df = stat_df[stat_df['Monotonic'] == True]

################### Trying more than two viral load groups ####################
all_par_ests = []
all_group_fits = []
x_vals = stat_df['Dist_X_Time'].describe()

#Get the dataframe for everyone except the current participant
low_group = stat_df[stat_df['Ave_VL'].lt(LOWER_THRESHOLD)]
low_group['Group'] = 'Low_VL'
low_group['Threshold'] = 0
mid_group = stat_df[stat_df['Ave_VL'].between(LOWER_THRESHOLD, UPPER_THRESHOLD)]
mid_group['Group'] = 'Mid_VL'
mid_group['Threshold'] = LOWER_THRESHOLD
high_group = stat_df[stat_df['Ave_VL'].gt(UPPER_THRESHOLD)]
high_group['Group'] = 'High_VL'
high_group['Threshold'] = UPPER_THRESHOLD

stat_df = pd.concat([low_group, mid_group, high_group], ignore_index=True)

#Estimate for the specific group
for name, group in stat_df.groupby('Group'):
    #Get the current estimate
    lower_fit, upper_fit, estimate_df = plne.bootstrap_rho(group,
                                                            NUM_BOOTSTRAPS)
    estimate_df['Threshold'] = group['Threshold'].unique()[0]
    estimate_df['Group'] = name

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
    group_fits['Threshold'] =  group['Threshold'].unique()[0]
    group_fits['Group'] = name

    all_group_fits.append(group_fits)
    all_par_ests.append(estimate_df)

all_par_ests = pd.concat(all_par_ests, ignore_index=True)
all_group_fits = pd.concat(all_group_fits, ignore_index=True)

all_conf_ints = []

#Calculate the confidence intervals for the estimates with low + high VL
for name, group in all_par_ests.groupby(['Group', 'Threshold']):
    lower_conf = np.quantile(group['Estimated_Rho'], 0.025)
    upper_conf = np.quantile(group['Estimated_Rho'], 0.975)

    #append the confidence interval
    all_conf_ints.append([lower_conf, upper_conf, name[0], name[1]])

all_conf_ints = pd.DataFrame(all_conf_ints, 
    columns=['lower_conf', 'upper_conf', 'Group', 'Threshold'])

############################# Plotting Panel 1 ################################
print(all_group_fits)
rcParams['mathtext.fontset'] = 'custom'
rcParams['mathtext.rm'] = 'DejaVu Sans'
rcParams['mathtext.it'] = 'DejaVu Sans:italic'
sns.set(rc={'figure.figsize':(20,7.5)}, font_scale = 2)
sns.set_palette("tab10")
fig_4, (ax1, ax2) = plt.subplots(1,2)

sns.boxplot(x ='Threshold', y ='Estimated_Rho', hue = 'Group', data = all_par_ests, ax = ax1)
ax1.set_xlabel(r'Group Viral Load Threshold (copies/ml)')
ax1.set_ylabel(r'Estimated Value of $\rho$')
ax1.ticklabel_format(style = 'sci', axis = 'y')
ax1.axhline(0.000008, linestyle = 'dashed', color = 'tab:green')
ax1.axhline(0.000014, color = 'tab:green')
ax1.axhline(0.00002, linestyle = 'dashed', color = 'tab:green')

# #plot the fits for the 50000 copies/mL threshold
low_50k = all_group_fits[all_group_fits['Group'] == 'Low_VL']
mid_50k = all_group_fits[all_group_fits['Group'] == 'Mid_VL']
high_50k = all_group_fits[all_group_fits['Group'] == 'High_VL']
sns.lineplot(x ='bin_edges', y ='ratio_bins', data = all_group_fits, hue = 'Group', ax = ax2)
sns.lineplot(x ='bin_edges', y ='mid_conf', data = all_group_fits, hue = 'Group', linestyle = 'dashed', ax = ax2)
sns.lineplot(x ='bin_edges', y ='lower_conf', data = all_group_fits, hue = 'Group', linestyle = 'dashed', ax = ax2)
sns.lineplot(x ='bin_edges', y ='high_conf', data = all_group_fits, hue = 'Group', linestyle = 'dashed', ax = ax2)

ax2.set_xlabel(r'Distance X Time (bp/generation)')
ax2.set_ylabel("D\' Ratio")


plt.tight_layout()
plt.savefig(outDir + 'threegroups.jpg')
plt.close()
