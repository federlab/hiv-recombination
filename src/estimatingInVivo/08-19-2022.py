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

NUM_BOOTSTRAPS = 1000

#For running on desktop
dataDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/"
vlDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/"
outDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/zanini/08-19-2022/"

#Today I am going to try comparing the data for fragment 4 and fragment 6
#I will start by plotting the viral load data for each fragment

#Make the dataframe containg D' ratios
stat_df = zu.combine_drats(dataDir)
stat_df = zu.label_vl_drats(stat_df, vlDir)

#Separate the dataframes for each fragment
stat_df_f4 = stat_df[stat_df['Fragment'] == 'F4']
stat_df_f6 = stat_df[stat_df['Fragment'] == 'F6']

####################### Plot the viral load data ##############################
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5), sharex= True, sharey=True)
sns.histplot(data=stat_df_f4, x='Ave_VL', ax = ax1, bins = 100)
ax1.set_xlabel('Viral Load F4 (copies/mL)')
sns.histplot(data=stat_df_f6, x='Ave_VL', ax = ax2, bins = 100)
ax2.set_xlabel('Viral Load F6 (copies/mL)')
plt.savefig(outDir + 'vl_hist.png')
plt.close()

####################### Get the Estimate for fragment 4 #######################
all_par_ests = []
all_group_fits = []
all_par_ests = []

#Get the estimate for fragment 4
lower_fit, upper_fit, estimate_df = plne.bootstrap_rho(stat_df_f4,
                                                        NUM_BOOTSTRAPS)
estimate_df['Fragment'] = 'F4'

#Bin the d' ratios so they are easier to view on the plots
binned_rat, binedges, bin_nums = stats.binned_statistic(
    stat_df_f4['Dist_X_Time'].to_numpy(), 
    stat_df_f4['d_ratio'].to_numpy(), bins = 100)

#Get the mean value of the estimates
mid_fit = np.quantile(estimate_df['Estimated_Rho'], 0.5)
mid_fit = estimate_df.iloc[(estimate_df['Estimated_Rho'] - mid_fit).abs().argsort()[:2]]

#Get the low and high values for the confidence interval
fit_vals_low = [plne.neher_leitner(x, lower_fit['c0'].to_numpy()[0], lower_fit['c1'].to_numpy()[0], lower_fit['c2'].to_numpy()[0]) for x in binedges]
fit_vals_high = [plne.neher_leitner(x, upper_fit['c0'].to_numpy()[0], upper_fit['c1'].to_numpy()[0], upper_fit['c2'].to_numpy()[0]) for x in binedges]
fit_vals_mid = [plne.neher_leitner(x, mid_fit['c0'].to_numpy()[0], mid_fit['c1'].to_numpy()[0], mid_fit['c2'].to_numpy()[0]) for x in binedges]

#Gather all of the fits for the participant
group_fits = pd.DataFrame(zip(binned_rat, binedges, fit_vals_low, fit_vals_mid, fit_vals_high),
                columns=['ratio_bins', 'bin_edges', 'lower_conf', 'mid_conf', 'high_conf'])
group_fits['Fragment'] = 'F4'

all_group_fits.append(group_fits)
all_par_ests.append(estimate_df)

####################### Get the Estimate for fragment 6 #######################
#Get the estimate for fragment 4
lower_fit, upper_fit, estimate_df = plne.bootstrap_rho(stat_df_f6,
                                                        NUM_BOOTSTRAPS)
estimate_df['Fragment'] = 'F6'

#Bin the d' ratios so they are easier to view on the plots
binned_rat, binedges, bin_nums = stats.binned_statistic(
    stat_df_f6['Dist_X_Time'].to_numpy(), 
    stat_df_f6['d_ratio'].to_numpy(), bins = 100)

#Get the mean value of the estimates
mid_fit = np.quantile(estimate_df['Estimated_Rho'], 0.5)
mid_fit = estimate_df.iloc[(estimate_df['Estimated_Rho'] - mid_fit).abs().argsort()[:2]]

#Get the low and high values for the confidence interval
fit_vals_low = [plne.neher_leitner(x, lower_fit['c0'].to_numpy()[0], lower_fit['c1'].to_numpy()[0], lower_fit['c2'].to_numpy()[0]) for x in binedges]
fit_vals_high = [plne.neher_leitner(x, upper_fit['c0'].to_numpy()[0], upper_fit['c1'].to_numpy()[0], upper_fit['c2'].to_numpy()[0]) for x in binedges]
fit_vals_mid = [plne.neher_leitner(x, mid_fit['c0'].to_numpy()[0], mid_fit['c1'].to_numpy()[0], mid_fit['c2'].to_numpy()[0]) for x in binedges]

#Gather all of the fits for the participant
group_fits = pd.DataFrame(zip(binned_rat, binedges, fit_vals_low, fit_vals_mid, fit_vals_high),
                columns=['ratio_bins', 'bin_edges', 'lower_conf', 'mid_conf', 'high_conf'])
group_fits['Fragment'] = 'F6'

all_group_fits.append(group_fits)
all_par_ests.append(estimate_df)

################# Put all of the data into the same dataframe #################
all_par_ests = pd.concat(all_par_ests, ignore_index=True)
all_group_fits = pd.concat(all_group_fits, ignore_index=True)

all_conf_ints = []

#Calculate the confidence intervals for the estimates with low + high VL
for name, group in all_par_ests.groupby(['Fragment']):
    lower_conf = np.quantile(group['Estimated_Rho'], 0.025)
    upper_conf = np.quantile(group['Estimated_Rho'], 0.975)

    #append the confidence interval
    all_conf_ints.append([lower_conf, upper_conf, name])

all_conf_ints = pd.DataFrame(all_conf_ints, 
    columns=['lower_conf', 'upper_conf', 'Fragment'])

print(all_par_ests)
print(all_conf_ints)

#Plot the fits for each fragment
data_F4 = all_group_fits[all_group_fits['Fragment'] == 'F4']
data_F6 = all_group_fits[all_group_fits['Fragment'] == 'F6']

sns.lineplot(x ='bin_edges', y ='ratio_bins', data = data_F4, color = 'tab:blue', label = 'F4')
sns.lineplot(x ='bin_edges', y ='mid_conf', data = data_F4, color = 'tab:blue', linestyle = 'dashed')
sns.lineplot(x ='bin_edges', y ='lower_conf', data = data_F4, color = 'navy', linestyle = 'dashed')
sns.lineplot(x ='bin_edges', y ='high_conf', data = data_F4, color = 'navy', linestyle = 'dashed')

sns.lineplot(x ='bin_edges', y ='ratio_bins', data = data_F6, color = 'tab:orange', label = 'F6')
sns.lineplot(x ='bin_edges', y ='mid_conf', data = data_F6, color = 'tab:orange', linestyle = 'dashed')
sns.lineplot(x ='bin_edges', y ='lower_conf', data = data_F6, color = 'saddlebrown', linestyle = 'dashed')
sns.lineplot(x ='bin_edges', y ='high_conf', data = data_F6, color = 'saddlebrown', linestyle = 'dashed')
ax2.set_xlabel(r'Distance X Time (bp/generation)')
ax2.set_ylabel("D\' Ratio")
plt.legend()
plt.savefig(outDir + "Fit_Results.png")
plt.close()