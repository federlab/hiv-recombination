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

GROUP_THRESHOLD_LIST = [7500, 10000, 25000, 50000, 100000, 200000]
NUM_BOOTSTRAPS = 1000
PAR_LIST = ['p9']


#For running on desktop
dataDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/"
vlDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/"
outDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/zanini/08-01-2022/"

#Make the dataframe containg D' ratios
stat_df = zu.combine_drats(dataDir)
stat_df = zu.label_vl_drats(stat_df, vlDir)

stat_df = stat_df[stat_df['Fragment'] != 'F5']
stat_df = stat_df[stat_df['Time_Diff'].gt(50)]

#Only get values for which the viral load doesn't leave the boundaries of the 
#initial and final measurments
stat_df = stat_df[stat_df['Monotonic'] == True]

############################# Panel 1 Analysis ################################
all_par_ests = []
all_group_fits = []
x_vals = stat_df['Dist_X_Time'].describe()

#Estimate rates specifically excluding each individual
for curr_thresh in GROUP_THRESHOLD_LIST:
    print(curr_thresh)
    for curr_par in PAR_LIST:
        print(curr_par)
        #Get the dataframe for everyone except the current participant
        stat_df['High_VL'] = stat_df['Ave_VL'].gt(curr_thresh)
        curr_stat_df = stat_df[stat_df['Participant'] == curr_par]

        #Estimate for the specific group
        for name, group in curr_stat_df.groupby('High_VL'):
            #Name the group by whether it is high or low viral load
            if name:
                curr_name = 'High_VL'
            else: curr_name = 'Low_VL'

            #Get the current estimate
            lower_fit, upper_fit, estimate_df = plne.bootstrap_rho(group,
                                                                NUM_BOOTSTRAPS)
            estimate_df['Threshold'] = curr_thresh
            estimate_df['Group'] = curr_name
            estimate_df['Participant'] = curr_par

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
            group_fits['Participant'] = curr_par

            all_group_fits.append(group_fits)
            all_par_ests.append(estimate_df)

all_par_ests = pd.concat(all_par_ests, ignore_index=True)
all_group_fits = pd.concat(all_group_fits, ignore_index=True)

all_conf_ints = []

#Calculate the confidence intervals for the estimates with low + high VL
for name, group in all_par_ests.groupby(['Group', 'Threshold', 'Participant']):
    lower_conf = np.quantile(group['Estimated_Rho'], 0.025)
    upper_conf = np.quantile(group['Estimated_Rho'], 0.975)

    #append the confidence interval
    all_conf_ints.append([lower_conf, upper_conf, name[0], name[1], name[2]])

all_conf_ints = pd.DataFrame(all_conf_ints, 
    columns=['lower_conf', 'upper_conf', 'Group', 'Threshold', 'Participant'])

print(all_conf_ints)

############################# Plotting Panel 1 ################################
#plot the 

# print(all_group_fits.columns)
# print(all_par_ests.columns)
rcParams['mathtext.fontset'] = 'custom'
rcParams['mathtext.rm'] = 'DejaVu Sans'
# rcParams['mathtext.it'] = 'DejaVu Sans:italic'
sns.set(font_scale = 2)
sns.set_palette("tab10")


myplot = sns.FacetGrid(all_par_ests, col = 'Participant', height = 10, aspect = 1.5)
myplot.map_dataframe(sns.boxplot, x ='Threshold', y ='Estimated_Rho', hue = 'Group')
plt.ylim(0, 0.0003)
plt.xlabel(r'Group Viral Load Threshold (copies/ml)')
plt.ylabel(r'Estimated Value of $\rho$')
plt.ticklabel_format(style = 'sci', axis = 'y')
plt.axhline(0.000008, linestyle = 'dashed', color = 'tab:green')
plt.axhline(0.000014, color = 'tab:green')
plt.axhline(0.00002, linestyle = 'dashed', color = 'tab:green')

# # #plot the fits for the 50000 copies/mL threshold
# data_50k = all_group_fits[all_group_fits['Threshold'] == 50000]
# low_50k = data_50k[data_50k['Group'] == 'Low_VL']
# high_50k = data_50k[data_50k['Group'] == 'High_VL']
# sns.lineplot(x ='bin_edges', y ='ratio_bins', data = low_50k, color = 'tab:blue', ax = ax2)
# sns.lineplot(x ='bin_edges', y ='mid_conf', data = low_50k, color = 'tab:blue', linestyle = 'dashed', ax = ax2)
# sns.lineplot(x ='bin_edges', y ='lower_conf', data = low_50k, color = 'gray', linestyle = 'dashed', ax = ax2)
# sns.lineplot(x ='bin_edges', y ='high_conf', data = low_50k, color = 'gray', linestyle = 'dashed', ax = ax2)

# sns.lineplot(x ='bin_edges', y ='ratio_bins', data = high_50k, color = 'tab:orange', ax = ax2)
# sns.lineplot(x ='bin_edges', y ='mid_conf', data = high_50k, color = 'tab:orange', linestyle = 'dashed', ax = ax2)
# sns.lineplot(x ='bin_edges', y ='lower_conf', data = high_50k, color = 'black', linestyle = 'dashed', ax = ax2)
# sns.lineplot(x ='bin_edges', y ='high_conf', data = high_50k, color = 'black', linestyle = 'dashed', ax = ax2)
# ax2.set_xlabel(r'Distance X Time (bp/generation)')
# ax2.set_ylabel("D\' Ratio")

plt.legend()
plt.tight_layout()
plt.savefig(outDir + 'within_host_results_1.jpg')
plt.close()

