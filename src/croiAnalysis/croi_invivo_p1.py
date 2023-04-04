import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
#for running on desktop
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plot_neher as plne
import zaniniUtil as zu
from scipy import stats
from matplotlib import rcParams

GROUP_THRESHOLD_LIST = [10000, 25000, 50000, 100000, 200000]
NUM_BOOTSTRAPS = 1000
EXAMPLE_THRESHOLD = 50000
PARTICIPANT = 'p1'


#For running on desktop
dataDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/"
vlDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/"
outDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/croiAnalysis/"

# #For running on cluster
# dataDir = "/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/"
# vlDir = "/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/"
# outDir = "/net/feder/vol1/home/evromero/2021_hiv-rec/results/croiAnalysis/"

#Make the dataframe containg D' ratios
stat_df = zu.combine_drats(dataDir)
stat_df = zu.label_vl_drats(stat_df, vlDir)

stat_df = stat_df[stat_df['Fragment'] != 'F5']
stat_df = stat_df[stat_df['Time_Diff'].gt(50)]
stat_df = stat_df[stat_df['Participant'] == PARTICIPANT]

#Only get values for which the viral load doesn't leave the boundaries of the 
#initial and final measurments
stat_df = stat_df[stat_df['Monotonic'] == True]

stat_df = stat_df[stat_df['d_i'] > 0.2]

############################# Panel 1 Analysis ################################
all_par_ests = []
all_group_fits = []
group_size_df = []
x_vals = stat_df['Dist_X_Time'].describe()


for curr_thresh in GROUP_THRESHOLD_LIST:

    #Get the dataframe for only the current participant
    stat_df['High_VL'] = stat_df['Ave_VL'].gt(curr_thresh)

    #Estimate for the specific group
    for name, group in stat_df.groupby('High_VL'):
        #Name the group by whether it is high or low viral load
        if name:
            curr_name = 'High_VL'
        else: curr_name = 'Low_VL'

        #Get the current estimate
        lower_fit, upper_fit, estimate_df = plne.bootstrap_rho(group,
                                                            NUM_BOOTSTRAPS)
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

        #Add the size of the group to the dictionary which labels the axis
        group_size_df.append([curr_thresh, curr_name, group.shape[0]])

        all_group_fits.append(group_fits)
        all_par_ests.append(estimate_df)

all_par_ests = pd.concat(all_par_ests, ignore_index=True)
all_group_fits = pd.concat(all_group_fits, ignore_index=True)
all_group_sizes = pd.DataFrame(group_size_df, columns=['Threshold', 'Group', 'Size'])

all_conf_ints = []

#Calculate the confidence intervals for the estimates with low + high VL
for name, group in all_par_ests.groupby(['Group', 'Threshold']):
    lower_conf = np.quantile(group['Estimated_Rho'], 0.025)
    mid_conf = np.quantile(group['Estimated_Rho'], 0.5)
    upper_conf = np.quantile(group['Estimated_Rho'], 0.975)

    #append the confidence interval
    all_conf_ints.append([lower_conf, mid_conf, upper_conf, name[0], name[1]])

all_conf_ints = pd.DataFrame(all_conf_ints, 
    columns=['lower_conf', 'mid_conf', 'upper_conf', 'Group', 'Threshold'])


#Output the data that made the figures
all_par_ests.to_csv(outDir + 'all_par_ests_p1_'+ str(NUM_BOOTSTRAPS)+ ' .csv')
all_conf_ints.to_csv(outDir + 'all_conf_ints_p1_' + str(NUM_BOOTSTRAPS)+ '.csv')
print(all_conf_ints)

############################# Plotting the Estimates ################################


print(all_group_fits)
rcParams.update({'font.size': 20,'figure.figsize':(15, 7.5) })
linewidth = 2

fig, (ax1, ax2) = plt.subplots(1, 2)
plt.subplots_adjust(hspace=0.25)

curr_data = all_par_ests[all_par_ests['Threshold'] == EXAMPLE_THRESHOLD]
sns.boxplot(x ='Threshold', y ='Estimated_Rho', hue = 'Group', data = curr_data, ax = ax2,  fliersize = 2, linewidth = linewidth)

ax2.set_ylabel(r'Estimated Recombination Rate ($\hat{\rho}$)')
ax2.set_xlabel(r'Group Viral Load Threshold (copies/ml)')
ax2.set_ylim(0, 0.00025)
ax2.ticklabel_format(style = 'sci', axis = 'y', scilimits = (0,0))
ax2.axhline(0.000008, linestyle = 'dashed', color = 'tab:green', linewidth = linewidth)
ax2.axhline(0.000014, color = 'tab:green', linewidth = linewidth)
ax2.axhline(0.00002, linestyle = 'dashed', color = 'tab:green', linewidth = linewidth)
ax2.xaxis.set_tick_params(labelbottom=True)
new_labels = ['Low VL', 'High VL']
for t, l in zip(ax2.get_legend().texts, new_labels):
    t.set_text(l)

####################### Plotting The Viral Loads######################################

available_files = os.listdir(vlDir)

#data frame of viral loads for each patient
viralLoadData = zu.make_viral_load_df(vlDir)
viralLoadData = viralLoadData[viralLoadData['Participant'] == 'p1']
viralLoadData['BelowThreshold'] = ['Low_VL' if x < EXAMPLE_THRESHOLD else 'High_VL' for x in viralLoadData['Viral load [virions/ml]']]

#now plot the viral loads
rcParams.update({'font.size': 18, 'figure.figsize':(13, 7)})

myplot = sns.lineplot(x = 'Days from infection', y = 'Viral load [virions/ml]', hue = 'BelowThreshold', style = 'BelowThreshold',data = viralLoadData,
                     errorbar = None, ax = ax1, markers = True, linewidth = linewidth)

ax1.axhline(EXAMPLE_THRESHOLD, linestyle = 'dashed', color = 'black', linewidth = linewidth)
ax1.set_yscale('log')
ax1.set_xlabel("Days from Infection")
ax1.set_ylabel("Viral Load [virions/ml]")
ax1.legend(loc = 'upper left')
plt.tight_layout()
plt.savefig(outDir + 'within_ests_invivo_' + str(NUM_BOOTSTRAPS) + '.jpg', dpi = 300)
plt.close()


####################### Plotting The Fit######################################
rcParams.update({'font.size': 18,'figure.figsize':(12, 5) })
# #plot the fits for the 50000 copies/mL threshold
fig, axs = plt.subplots(1, 2)

data_50k = all_group_fits[all_group_fits['Threshold'] == EXAMPLE_THRESHOLD]
low_50k = data_50k[data_50k['Group'] == 'Low_VL']
high_50k = data_50k[data_50k['Group'] == 'High_VL']
sns.lineplot(x ='bin_edges', y ='ratio_bins', data = low_50k, color = 'tab:blue', linewidth = linewidth, label = 'Low VL', ax = axs[1])
sns.lineplot(x ='bin_edges', y ='mid_conf', data = low_50k, color = 'navy',linewidth = linewidth, ax = axs[1])
sns.lineplot(x ='bin_edges', y ='lower_conf', data = low_50k, color = 'navy', linestyle = 'dashed', linewidth = linewidth, ax = axs[1])
sns.lineplot(x ='bin_edges', y ='high_conf', data = low_50k, color = 'navy', linestyle = 'dashed', linewidth = linewidth, ax = axs[1])

sns.lineplot(x ='bin_edges', y ='ratio_bins', data = high_50k, color = 'tab:orange', linewidth = linewidth, label = 'High VL', ax = axs[1])
sns.lineplot(x ='bin_edges', y ='mid_conf', data = high_50k, color = 'saddlebrown', linewidth = linewidth, ax = axs[1])
sns.lineplot(x ='bin_edges', y ='lower_conf', data = high_50k, color = 'saddlebrown', linestyle = 'dashed', linewidth = linewidth, ax = axs[1])
sns.lineplot(x ='bin_edges', y ='high_conf', data = high_50k, color = 'saddlebrown', linestyle = 'dashed', linewidth = linewidth, ax = axs[1])
axs[1].get_legend().remove()
axs[1].set_xlabel(r'Distance $\cdot$ Time (bp $\cdot$ generation)')
axs[1].set_ylabel("D\' Ratio")


####################### Plotting The Viral Loads######################################

available_files = os.listdir(vlDir)

#data frame of viral loads for each patient
viralLoadData = zu.make_viral_load_df(vlDir)
viralLoadData = viralLoadData[viralLoadData['Participant'] == 'p1']
viralLoadData['BelowThreshold'] = ['Low_VL' if x < EXAMPLE_THRESHOLD else 'High_VL' for x in viralLoadData['Viral load [virions/ml]']]

#now plot the viral loads
rcParams.update({'font.size': 18, 'figure.figsize':(13, 7)})

myplot = sns.lineplot(x = 'Days from infection', y = 'Viral load [virions/ml]', hue = 'BelowThreshold', style = 'BelowThreshold',data = viralLoadData,
                     errorbar = None, ax = axs[0], markers = True)

axs[0].axhline(EXAMPLE_THRESHOLD, linestyle = 'dashed', color = 'black', linewidth = linewidth)
axs[0].set_yscale('log')
axs[0].set_xlabel("Days from Infection")
axs[0].set_ylabel("Viral Load [virions/ml]")
axs[0].legend(loc = 'upper left')
plt.tight_layout()
plt.savefig(outDir + 'within_p1_fits_invivo_' + str(NUM_BOOTSTRAPS) + '.jpg', dpi = 300)
plt.close()