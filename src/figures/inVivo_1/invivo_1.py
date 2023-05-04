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
QUANTILE_LIST = [0, 0.33, 0.66, 1]

#For running on desktop
dataDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/"
vlDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/"
outDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/paper/inVivo_1/"

#For running on cluster
# dataDir = "/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/"
# vlDir = "/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/"
# outDir = "/net/feder/vol1/home/evromero/2021_hiv-rec/results/paper/inVivo_1/"

#Make the dataframe containing D' ratios
stat_df = zu.combine_drats(dataDir)
stat_df = zu.label_vl_drats(stat_df, vlDir)

stat_df = stat_df[stat_df['Fragment'] != 'F5']
stat_df = stat_df[stat_df['Time_Diff'].gt(50)]

#Only get values for which the viral load doesn't leave the boundaries of the 
#initial and final measurments
stat_df = stat_df[stat_df['Monotonic'] == True]
print(min(stat_df['d_i']))

######################### Configure Plot Settings #############################
#plot the estimates to show how accurate they are
params = {'figure.figsize':(7, 2.5), 'axes.labelsize': 6,'axes.titlesize':6,  
          'legend.fontsize': 6, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
          'legend.title_fontsize': 6}
rcParams.update(params)

#Make a three panel figure
fig, axs = plt.subplots(1, 3)
linewidth = 1
markersize = 4


###################### Analysis for panels 2 and 3 ############################
all_par_ests = []
all_group_fits = []
group_size_df = []

lower_label = "Lower \n<" + str(int(np.quantile(stat_df['Ave_VL'], 0.33)))
middle_label = "Middle \n" + str(int(np.quantile(stat_df['Ave_VL'], 0.33))) + " - " + str(int(np.quantile(stat_df['Ave_VL'], 0.66))) 
upper_label = "Upper \n>" + str(int(np.quantile(stat_df['Ave_VL'], 0.66)))

#First we need to group the data by the viral load quantiles
stat_df['VL_Quantile'], bins = pd.qcut(stat_df['Ave_VL'], q = QUANTILE_LIST, retbins= True, labels = [lower_label, middle_label, upper_label])
grouped_stats = stat_df.groupby(by = ['VL_Quantile'])
print(bins)



for name, group in grouped_stats:
    #Get the current estimate
    lower_fit, upper_fit, estimate_df = plne.bootstrap_rho(group,
                                                            NUM_BOOTSTRAPS)
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
    group_fits['Group'] = name

    #Add the size of the group to the dictionary which labels the axis
    group_size_df.append([name, group.shape[0]])
    
    #Add all of the current results
    all_par_ests.append(estimate_df)
    all_group_fits.append(group_fits)

all_par_ests = pd.concat(all_par_ests, ignore_index=True)
all_group_fits = pd.concat(all_group_fits, ignore_index=True)
all_group_sizes = pd.DataFrame(group_size_df, columns=['Group', 'Size'])
print(all_group_sizes)
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
all_par_ests.to_csv(outDir + 'all_par_ests' + str(NUM_BOOTSTRAPS) + '.csv')
all_conf_ints.to_csv(outDir + 'all_conf_ints' + str(NUM_BOOTSTRAPS) + '.csv')

###################### Panel 1: Viral Load vs. Time ###########################
#data frame of viral loads for each patient
viralLoadData = zu.make_viral_load_df(vlDir)

#now plot the viral loads for each participant
par_list = viralLoadData['Participant'].unique()
for curr_par in par_list:
    curr_vlData = viralLoadData.copy(deep = True)
    curr_vlData = curr_vlData[curr_vlData['Participant'] == curr_par]
    # curr_vlData['VL Group'] = pd.cut(curr_vlData['Viral load [virions/ml]'], bins)

    myplot = sns.lineplot(x = 'Days from infection', 
                          y = 'Viral load [virions/ml]', data = curr_vlData,
                        color = 'k', errorbar = None, ax = axs[0], markers = True, 
                        linewidth = linewidth, markersize = markersize)

#plot the quantile thresholds
color_list = ['tab:blue', 'tab:gray', 'tab:orange']

aug_bins = [1] + bins.tolist()[1:-1] + [10**6]
fill_x = np.linspace(0, 6000, 100)

for i in range(0, len(QUANTILE_LIST)-1):
    axs[0].fill_between(fill_x, aug_bins[i], aug_bins[i+1], color = color_list[i], alpha = 0.3)
axs[0].set_yscale('log')
axs[0].set_ylim(np.min(viralLoadData['Viral load [virions/ml]']), np.max(viralLoadData['Viral load [virions/ml]']))
axs[0].set_xlim(0, 6000)
axs[0].set_xlabel("Time (days since EDI)")
axs[0].set_ylabel("Viral Load (copies/mL)")


###################### Panel 2: D' Ratio vs. dDeltaT ##########################

#Plot the lineplots 
sns.lineplot(x ='bin_edges', y ='ratio_bins', data = all_group_fits, hue = 'Group', ax = axs[1], linewidth = linewidth,
             palette=['tab:blue', 'tab:gray', 'tab:orange'])
sns.lineplot(x ='bin_edges', y ='mid_conf', data = all_group_fits, hue = 'Group', ax = axs[1], linewidth = linewidth,
             palette=['tab:blue', 'tab:gray', 'tab:orange'])

handles, labels = axs[1].get_legend_handles_labels()
labels = ['Lower', 'Middle', 'Upper']
by_label = dict(zip(labels, handles))
axs[1].legend(by_label.values(), by_label.keys(), title = 'Viral Load Tertile')


axs[1].set_xlabel(r'Distance $\cdot$ Time (bp $\cdot$ generations)')
axs[1].set_ylabel("D\' Ratio")


####################### Panel 3: Boxplot of Quartiles #########################

sns.boxplot(x ='Group', y ='Estimated_Rho', data = all_par_ests, ax = axs[2], fliersize = 2, linewidth = linewidth,
            palette=['tab:blue', 'tab:gray', 'tab:orange'])
axs[2].set_ylabel(r'Estimated Recombination Rate ($\hat{\rho}$)')
axs[2].set_xlabel(r'Viral Load Tertile (copies/mL)')

axs[2].ticklabel_format(style = 'sci', axis = 'y', scilimits = (0,0))
axs[2].xaxis.set_tick_params(labelbottom=True)

plt.tight_layout()
plt.savefig(outDir + 'invivo_1_' + str(NUM_BOOTSTRAPS) + '.jpg', dpi = 300)

plt.close()