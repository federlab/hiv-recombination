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

#For running on desktop
dataDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/"
vlDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/"
outDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/paper/fig3/"

#For running on cluster
# dataDir = "/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/"
# vlDir = "/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/"
# outDir = "/net/feder/vol1/home/evromero/2021_hiv-rec/results/paper/fig3/"

#Make the dataframe containing D' ratios
stat_df = zu.combine_drats(dataDir)
stat_df = zu.label_vl_drats(stat_df, vlDir)

stat_df = stat_df[stat_df['Fragment'] != 'F5']
stat_df = stat_df[stat_df['Time_Diff'].gt(50)]

#Only get values for which the viral load doesn't leave the boundaries of the 
#initial and final measurments
stat_df = stat_df[stat_df['Monotonic'] == True]
print(min(stat_df['d_i']))


############################# Panel 1 Analysis ################################
all_par_ests = []
all_group_fits = []
x_vals = stat_df['Dist_X_Time'].describe()
group_size_df = []

for curr_thresh in GROUP_THRESHOLD_LIST:
    #Get the dataframe for everyone except the current participant
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
        
        #Add all of the current results
        all_par_ests.append(estimate_df)
        all_group_fits.append(group_fits)


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
all_par_ests.to_csv(outDir + 'all_par_ests' + str(NUM_BOOTSTRAPS) + '.csv')
all_conf_ints.to_csv(outDir + 'all_conf_ints' + str(NUM_BOOTSTRAPS) + '.csv')
######################### Making the figure layout ############################

rcParams.update({'font.size': 8,})
linewidth = 1

fig, axd = plt.subplot_mosaic([['upper left', 'right'],
                               ['lower left', 'right'],
                               ['lower left', 'bottom right']], height_ratios=[0.75, 0.5, 0.25],
                              figsize=(6.5, 4), layout="constrained")



############################# Plotting Panel C ################################

sns.boxplot(x ='Threshold', y ='Estimated_Rho', hue = 'Group', data = all_par_ests, ax = axd['right'], palette = ['tab:blue', 'tab:orange'], fliersize = 2, linewidth = linewidth)
# print(axd['right'])
axd['right'].set_ylabel(r'Estimated Recombination Rate ($\hat{\rho}$)')
axd['right'].set_xlabel(r'Group Viral Load Threshold (copies/mL)')
axd['right'].ticklabel_format(style = 'sci', axis = 'y', scilimits = (0,0))
axd['right'].axhline(0.000008, linestyle = 'dashed', color = 'tab:green', linewidth = linewidth)
axd['right'].axhline(0.000014, color = 'tab:green', linewidth = linewidth)
axd['right'].axhline(0.00002, linestyle = 'dashed', color = 'tab:green', linewidth = linewidth)
axd['right'].set_ylim(0, 0.0003)
axd['right'].xaxis.set_tick_params(labelbottom=True)


# axd['right'].legend(bbox_to_anchor=(1, 0.5))
new_labels = ['Low VL', 'High VL']
for t, l in zip(axd['right'].get_legend().texts, new_labels):
    t.set_text(l)


############################# Plotting Panel D ################################
print(all_group_sizes)

sns.barplot(x = 'Threshold', y = 'Size', hue = 'Group', data = all_group_sizes, ci = None, ax = axd['bottom right'], palette = ['tab:blue', 'tab:orange'])
axd['bottom right'].set_xlabel(r'Group Viral Load Threshold (copies/mL)')
axd['bottom right'].set_ylabel("Group Size")
axd['bottom right'].get_legend().remove()
axd['bottom right'].ticklabel_format(style = 'sci', axis = 'y', scilimits = (0,0))
axd['bottom right'].xaxis.set_tick_params(labelbottom=True)
#fig.align_ylabels([axd['right'], axd['bottom right']])
plt.tight_layout()



####################### Plotting Panel B ######################################

# #plot the fits for the 50000 copies/mL threshold
data_50k = all_group_fits[all_group_fits['Threshold'] == EXAMPLE_THRESHOLD]
low_50k = data_50k[data_50k['Group'] == 'Low_VL']
high_50k = data_50k[data_50k['Group'] == 'High_VL']
sns.lineplot(x ='bin_edges', y ='ratio_bins', data = low_50k, color = 'tab:blue', ax = axd['lower left'], linewidth = linewidth)
sns.lineplot(x ='bin_edges', y ='mid_conf', data = low_50k, color = 'navy', ax = axd['lower left'],linewidth = linewidth)
sns.lineplot(x ='bin_edges', y ='lower_conf', data = low_50k, color = 'navy', linestyle = 'dashed', ax = axd['lower left'], linewidth = linewidth)
sns.lineplot(x ='bin_edges', y ='high_conf', data = low_50k, color = 'navy', linestyle = 'dashed', ax = axd['lower left'], linewidth = linewidth)

sns.lineplot(x ='bin_edges', y ='ratio_bins', data = high_50k, color = 'tab:orange', ax = axd['lower left'], linewidth = linewidth)
sns.lineplot(x ='bin_edges', y ='mid_conf', data = high_50k, color = 'saddlebrown', ax = axd['lower left'], linewidth = linewidth)
sns.lineplot(x ='bin_edges', y ='lower_conf', data = high_50k, color = 'saddlebrown', linestyle = 'dashed', ax = axd['lower left'], linewidth = linewidth)
sns.lineplot(x ='bin_edges', y ='high_conf', data = high_50k, color = 'saddlebrown', linestyle = 'dashed', ax = axd['lower left'], linewidth = linewidth)
axd['lower left'].set_xlabel(r'Distance $\cdot$ Time (bp $\cdot$ generations)')
axd['lower left'].set_ylabel("D\' Ratio")



############################# Plotting Panel A ################################
#now loop through pairs of timepoints and plot them on the grid
grouped_par = stat_df.groupby(by = ['Day_1', 'Day_2'])
for name, group in grouped_par:

    time_1 = name[0]
    time_2 = name[1]

    vl_1 = group['VL_1'].unique()[0]
    vl_2 = group['VL_2'].unique()[0]

    ave_vl = group['Ave_VL'].unique()[0]
    ave_time = np.mean([time_1, time_2])
    if  ave_vl > EXAMPLE_THRESHOLD:
        my_color = "tab:orange"
        my_color2 = "peru"
    else: 
        my_color = "tab:blue"
        my_color2 = "navy"

    axd['upper left'].plot([time_1, time_2], [vl_1, vl_2], color = my_color, linewidth = linewidth)
    axd['upper left'].axhline(EXAMPLE_THRESHOLD, color = "gray", linestyle = "--", linewidth = linewidth)

for name, group in grouped_par:

    time_1 = name[0]
    time_2 = name[1]

    vl_1 = group['VL_1'].unique()[0]
    vl_2 = group['VL_2'].unique()[0]

    ave_vl = group['Ave_VL'].unique()[0]
    ave_time = np.mean([time_1, time_2])
    if  ave_vl > EXAMPLE_THRESHOLD:
        my_color = "tab:orange"
        my_color2 = "saddlebrown"
    else: 
        my_color = "tab:blue"
        my_color2 = "navy"
    if time_1 > 4000:
        print(time_1, time_2, vl_1, vl_2, ave_vl, ave_time)


    axd['upper left'].plot(ave_time, ave_vl, color = my_color2, marker = "D", markersize = 2)

axd['upper left'].set_xlabel('Time (days since EDI)')
axd['upper left'].ticklabel_format(style = 'sci', axis = 'y', scilimits = (0,0))
axd['upper left'].set_ylabel('Viral Load (copies/mL)')
plt.tight_layout()
plt.savefig(outDir + 'fig3refit_' + str(NUM_BOOTSTRAPS) + '.jpg', dpi = 300)

plt.close()