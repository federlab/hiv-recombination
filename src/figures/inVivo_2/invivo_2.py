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

GROUP_THRESHOLD_LIST = [10000, 25000, 50000, 100000, 200000]
NUM_BOOTSTRAPS = 1000
PAR_LIST = ['p1']
EXAMPLE_THRESHOLD = 50000


#For running on desktop
dataDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/"
vlDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/"
outDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/paper/inVivo_2/"

# #For running on cluster
# dataDir = "/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/"
# vlDir = "/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/"
# outDir = "/net/feder/vol1/home/evromero/2021_hiv-rec/results/paper/inVivo_2/"

#Make the dataframe containg D' ratios
stat_df = zu.combine_drats(dataDir)
stat_df = zu.label_vl_drats(stat_df, vlDir)

stat_df = stat_df[stat_df['Fragment'] != 'F5']
stat_df = stat_df[stat_df['Time_Diff'].gt(50)]

#Only get values for which the viral load doesn't leave the boundaries of the 
#initial and final measurments
stat_df = stat_df[stat_df['Monotonic'] == True]

stat_df = stat_df[stat_df['d_i'] > 0.2]

############################# Panel 1 Analysis ################################
all_par_ests = []
all_group_fits = []
group_size_df = []
x_vals = stat_df['Dist_X_Time'].describe()


for curr_par in PAR_LIST:
    #Get the dataframe for only the current participant
    
    curr_stat_df = stat_df[stat_df['Participant'] == curr_par]
    curr_stat_df['VL_Quantile'], bins = pd.qcut(curr_stat_df['Ave_VL'], q = [0, 0.5, 1], retbins= True, labels = ["Low_VL","High_VL"])
    print(bins)
    sns.histplot(data = curr_stat_df, x = 'Ave_VL', hue = 'VL_Quantile', bins = 100)
    plt.savefig(outDir + 'fig4reFit_vlHist.jpg', dpi = 300)     
    #Estimate for the specific group
    for name, group in curr_stat_df.groupby('VL_Quantile'):
        #Name the group by whether it is high or low viral load
        curr_name = name

        #Get the current estimate
        lower_fit, upper_fit, estimate_df = plne.bootstrap_rho(group,
                                                            NUM_BOOTSTRAPS)
        estimate_df['Group'] = curr_name
        estimate_df['Participant'] = curr_par

        #Bin the d' ratios so they are easier to view on the plots
        binned_rat, binedges, bin_nums = stats.binned_statistic(
            group['Dist_X_Time'].to_numpy(), 
            group['d_ratio'].to_numpy(), bins = 50)

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
        group_fits['Participant'] = curr_par

        #Add the size of the group to the dictionary which labels the axis
        group_size_df.append([curr_name, group.shape[0]])

        all_group_fits.append(group_fits)
        all_par_ests.append(estimate_df)

all_par_ests = pd.concat(all_par_ests, ignore_index=True)
all_group_fits = pd.concat(all_group_fits, ignore_index=True)
all_group_sizes = pd.DataFrame(group_size_df, columns=['Group', 'Size'])

all_conf_ints = []

#Calculate the confidence intervals for the estimates with low + high VL
for name, group in all_par_ests.groupby(['Group']):
    lower_conf = np.quantile(group['Estimated_Rho'], 0.025)
    mid_conf = np.quantile(group['Estimated_Rho'], 0.5)
    upper_conf = np.quantile(group['Estimated_Rho'], 0.975)

    #append the confidence interval
    all_conf_ints.append([lower_conf, mid_conf, upper_conf, name[0]])

all_conf_ints = pd.DataFrame(all_conf_ints, 
    columns=['lower_conf', 'mid_conf', 'upper_conf', 'Group'])


#Output the data that made the figures
all_par_ests.to_csv(outDir + 'all_par_ests_'+ str(NUM_BOOTSTRAPS)+ ' .csv')
all_conf_ints.to_csv(outDir + 'all_conf_ints_' + str(NUM_BOOTSTRAPS)+ '.csv')
print(all_conf_ints)





############################# We want 3 panels ################################
#plot the estimates to show how accurate they are
params = {'figure.figsize':(7, 2.25), 'axes.labelsize': 6,'axes.titlesize':6,  
          'legend.fontsize': 6, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
          'legend.title_fontsize': 6}
rcParams.update(params)
fig, axd = plt.subplot_mosaic([['left', 'middle', 'right']], height_ratios=[1],
                              figsize=(7, 2.5), layout="constrained")


print(all_group_fits)
rcParams.update({'font.size': 8})
linewidth = 1
markersize = 4
markeredgewidth = 1
sns.set_palette("tab10")
sns.set_style("white")


############################ Plot viral loads #################################
ax0 = axd['left']

#data frame of viral loads for each patient
viralLoadData = zu.make_viral_load_df(vlDir)
viralLoadData = viralLoadData[viralLoadData['Participant'] == 'p1']
viralLoadData['BelowThreshold'] = ['Low_VL' if x < bins[1] else 'High_VL' for x in viralLoadData['Viral load [virions/ml]']]

#now plot the viral loads
myplot = sns.lineplot(x = 'Days from infection', y = 'Viral load [virions/ml]', hue = 'BelowThreshold', style = 'BelowThreshold',data = viralLoadData,
                     errorbar = None, ax = ax0, markers = True, linewidth = linewidth, markersize = markersize)
ax0.axhline(y = bins[1], color = 'black', linestyle = 'dashed', linewidth = linewidth)
ax0.set_yscale('log')
ax0.set_xlabel("Time (days since EDI)")
ax0.set_ylabel("Viral Load (copies/mL)")
ax0.legend(loc = "lower right",title = 'Viral Load Quantile')
new_labels = ['Below Median', 'Above Median']
for t, l in zip(ax0.get_legend().texts, new_labels):
    t.set_text(l)


lgnd = ax0.get_legend()
for handle in lgnd.legendHandles:
    handle.set_markersize(markersize)
    handle.set_linewidth(linewidth)

lgnd.legendHandles[1].set_marker("x")

########################## Plot the box plots ################################
ax1 = axd['right']

lower_label = "Below Median \n< " + str(int(bins[1]))
upper_label = "Above Median \n> " + str(int(bins[1]))


label_dict = {'Low_VL': lower_label, 'High_VL': upper_label}
all_par_ests['Group'] = [label_dict[x] for x in all_par_ests['Group']]

sns.boxplot(x ='Group', y ='Estimated_Rho', data = all_par_ests, ax = ax1,  fliersize = 2, linewidth = linewidth,
            palette = 'tab10')

ax1.set_ylabel(r'Estimated Recombination Rate ($\hat{\rho}$)')
ax1.set_xlabel(r'Viral Load Quantile (copies/mL)')
ax1.ticklabel_format(style = 'sci', axis = 'y', scilimits = (0,0))
ax1.xaxis.set_tick_params(labelbottom=True)
# new_labels = ['Low VL', 'High VL']
# for t, l in zip(ax1.get_legend().texts, new_labels):
#     t.set_text(l)

############################# Plotting Panel D ################################
print(all_group_sizes)

ax2 = axd['middle']

low_50k = all_group_fits[all_group_fits['Group'] == 'Low_VL']
high_50k = all_group_fits[all_group_fits['Group'] == 'High_VL']
sns.lineplot(x ='bin_edges', y ='ratio_bins', data = low_50k, color = 'tab:blue', linewidth = linewidth, label = 'Low VL', ax = ax2)
sns.lineplot(x ='bin_edges', y ='mid_conf', data = low_50k, color = 'navy',linewidth = linewidth, ax = ax2)
# sns.lineplot(x ='bin_edges', y ='lower_conf', data = low_50k, color = 'navy', linestyle = 'dashed', linewidth = linewidth, ax = ax2)
# sns.lineplot(x ='bin_edges', y ='high_conf', data = low_50k, color = 'navy', linestyle = 'dashed', linewidth = linewidth, ax = ax2)

sns.lineplot(x ='bin_edges', y ='ratio_bins', data = high_50k, color = 'tab:orange', linewidth = linewidth, label = 'High VL', ax = ax2)
sns.lineplot(x ='bin_edges', y ='mid_conf', data = high_50k, color = 'saddlebrown', linewidth = linewidth, ax = ax2)
# sns.lineplot(x ='bin_edges', y ='lower_conf', data = high_50k, color = 'saddlebrown', linestyle = 'dashed', linewidth = linewidth, ax = ax2)
# sns.lineplot(x ='bin_edges', y ='high_conf', data = high_50k, color = 'saddlebrown', linestyle = 'dashed', linewidth = linewidth, ax = ax2)
ax2.get_legend().remove()
ax2.set_xlabel(r'Distance $\cdot$ Time (bp $\cdot$ generations)')
ax2.set_ylabel("D\' Ratio")

plt.tight_layout()
plt.savefig(outDir + 'fig4reFit_' + str(NUM_BOOTSTRAPS) + '.jpg', dpi = 300)
plt.close()