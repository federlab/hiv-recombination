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
from matplotlib import gridspec

GROUP_THRESHOLD_LIST = [7500, 10000, 25000, 50000, 100000, 200000]
NUM_BOOTSTRAPS = 1000
EXAMPLE_THRESHOLD = 50000

#This file will give us grouped estimates across all participants (not
# separated by viral load)

#For running on desktop
dataDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/"
vlDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/"
outDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/paper/pooled_estimation/"

#For running on cluster
# dataDir = "/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/"
# vlDir = "/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/"
# outDir = "/net/feder/vol1/home/evromero/2021_hiv-rec/results/paper/pooled_estimation/"

#Make the dataframe containing D' ratios
stat_df = zu.combine_drats(dataDir)
stat_df = zu.label_vl_drats(stat_df, vlDir)

stat_df = stat_df[stat_df['Fragment'] != 'F5']
stat_df = stat_df[stat_df['Time_Diff'].gt(50)]

#Only get values for which the viral load doesn't leave the boundaries of the 
#initial and final measurments
stat_df = stat_df[stat_df['Monotonic'] == True]
print(min(stat_df['d_i']))


all_par_ests = []
all_group_fits = []
x_vals = stat_df['Dist_X_Time'].describe()
group_size_df = []


#Get the current estimate
lower_fit, upper_fit, estimate_df = plne.bootstrap_rho(stat_df,
                                                    NUM_BOOTSTRAPS)


#Bin the d' ratios so they are easier to view on the plots
binned_rat, binedges, bin_nums = stats.binned_statistic(
    stat_df['Dist_X_Time'].to_numpy(), 
    stat_df['d_ratio'].to_numpy(), bins = 100)

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


all_conf_ints = []


lower_conf = np.quantile(estimate_df['Estimated_Rho'], 0.025)
mid_conf = np.quantile(estimate_df['Estimated_Rho'], 0.5)
upper_conf = np.quantile(estimate_df['Estimated_Rho'], 0.975)

all_conf_ints = pd.DataFrame([[lower_conf, mid_conf, upper_conf]], 
    columns=['lower_conf', 'mid_conf', 'upper_conf'])


#Output the data that made the figures
print(all_conf_ints)
all_conf_ints.to_csv(outDir + 'all_conf_ints' + str(NUM_BOOTSTRAPS) + '.csv')

# #plot the fits for the 50000 copies/mL threshold

sns.lineplot(x ='bin_edges', y ='ratio_bins', data = group_fits, color = 'k')
sns.lineplot(x ='bin_edges', y ='mid_conf', data = group_fits, color = 'purple')
sns.lineplot(x ='bin_edges', y ='lower_conf', data = group_fits, color = 'red', linestyle = 'dashed')
sns.lineplot(x ='bin_edges', y ='high_conf', data = group_fits, color = 'blue', linestyle = 'dashed')


plt.xlabel(r'Distance $\cdot$ Time (bp $\cdot$ generation)')
plt.ylabel("D\' Ratio")
plt.savefig(outDir + 'all_fits' + str(NUM_BOOTSTRAPS) + '.png', dpi = 300)