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

GROUP_THRESHOLD_LIST = [10000, 25000, 50000, 100000]
NUM_BOOTSTRAPS = 1000
PAR_LIST = ['p3', 'p9']


#For running on desktop
dataDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/"
vlDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/"
outDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/paper/supp_within/"

# #For running on cluster
# dataDir = "/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/"
# vlDir = "/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/"
# outDir = "/net/feder/vol1/home/evromero/2021_hiv-rec/results/paper/supp_within/"

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

#Estimate rates specifically excluding each individual
for curr_thresh in GROUP_THRESHOLD_LIST:
    print(curr_thresh)
    for curr_par in PAR_LIST:
        print(curr_par)
        #Get the dataframe for only the current participant
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

            #Add the size of the group to the dictionary which labels the axis
            group_size_df.append([curr_thresh, curr_name, group.shape[0], curr_par])

            all_group_fits.append(group_fits)
            all_par_ests.append(estimate_df)

all_par_ests = pd.concat(all_par_ests, ignore_index=True)
all_group_fits = pd.concat(all_group_fits, ignore_index=True)
all_group_sizes = pd.DataFrame(group_size_df, columns=['Threshold', 'Group', 'Size', 'Participant'])

all_conf_ints = []

#Calculate the confidence intervals for the estimates with low + high VL
for name, group in all_par_ests.groupby(['Group', 'Threshold', 'Participant']):
    lower_conf = np.quantile(group['Estimated_Rho'], 0.025)
    mid_conf = np.quantile(group['Estimated_Rho'], 0.5)
    upper_conf = np.quantile(group['Estimated_Rho'], 0.975)

    #append the confidence interval
    all_conf_ints.append([lower_conf, mid_conf, upper_conf, name[0], name[1], name[2]])

all_conf_ints = pd.DataFrame(all_conf_ints, 
    columns=['lower_conf', 'mid_conf', 'upper_conf', 'Group', 'Threshold', 'Participant'])


#Output the data that made the figures
all_par_ests.to_csv(outDir + 'all_par_ests.csv')
all_conf_ints.to_csv(outDir + 'all_conf_ints.csv')
print(all_conf_ints)

############################# Plotting Panel C ################################


rcParams.update({'font.size': 8,'figure.figsize':(6.5, 4) })
linewidth = 1
sns.set_palette("tab10")
sns.set_style("white")
fig, axs = plt.subplots(2, 2, sharex=True, gridspec_kw={'height_ratios': [7, 1]})
hue_order = ['Low_VL', 'High_VL']

for i in range(len(PAR_LIST)):
    curr_par_ests = all_par_ests[all_par_ests['Participant'] == PAR_LIST[i]]
    curr_group_sizes = all_group_sizes[all_group_sizes['Participant'] == PAR_LIST[i]]
    plt.subplots_adjust(hspace=0.25)

    sns.boxplot(x ='Threshold', y ='Estimated_Rho', hue = 'Group', hue_order = hue_order, data = curr_par_ests, ax = axs[0][i],  fliersize = 2, linewidth = linewidth)

    axs[0][i].set_title('Participant ' + str(PAR_LIST[i][1]), fontsize = 8)
    axs[0][i].set_ylabel(r'Estimated Recombination Rate ($\hat{\rho}$)')
    axs[0][i].set_xlabel(r'Group Viral Load Threshold (copies/ml)')
    axs[0][i].set_ylim(0, 0.0003)
    axs[0][i].ticklabel_format(style = 'sci', axis = 'y', scilimits = (0,0))
    axs[0][i].axhline(0.000008, linestyle = 'dashed', color = 'tab:green', linewidth = linewidth)
    axs[0][i].axhline(0.000014, color = 'tab:green', linewidth = linewidth)
    axs[0][i].axhline(0.00002, linestyle = 'dashed', color = 'tab:green', linewidth = linewidth)
    axs[0][i].xaxis.set_tick_params(labelbottom=True)
    new_labels = ['Low VL', 'High VL']
    for t, l in zip(axs[0][i].get_legend().texts, new_labels):
        t.set_text(l)

    ############################# Plotting Panel D ################################
    #calculate the size of each group to help format bar labels
    print(curr_group_sizes)
    low_vl_labels = curr_group_sizes[curr_group_sizes['Group'] == 'Low_VL']
    low_vl_labels = low_vl_labels.copy()
    low_vl_labels.sort_values(by = 'Threshold', inplace = True)
    low_vl_labels = low_vl_labels['Size'].tolist()

    low_vl_labels = [x  if x < 1000 else '' for x in low_vl_labels]

 

    high_vl_labels = curr_group_sizes[curr_group_sizes['Group'] == 'High_VL']
    high_vl_labels = high_vl_labels.copy()
    high_vl_labels.sort_values(by = 'Threshold', inplace = True)
    high_vl_labels = high_vl_labels['Size'].tolist()

    high_vl_labels = [x  if x < 1000 else '' for x in high_vl_labels]
    print(low_vl_labels)
    print(high_vl_labels)
    


    sns.barplot(x = 'Threshold', y = 'Size', hue = 'Group', hue_order = hue_order, data = curr_group_sizes, errorbar = None, ax = axs[1][i])
    if i == 1:
        axs[1][i].bar_label(axs[1][i].containers[0], labels = low_vl_labels, fontsize = 5)
        axs[1][i].bar_label(axs[1][i].containers[1], labels = high_vl_labels, fontsize = 5)
    axs[1][i].set_xlabel(r'Group Viral Load Threshold (copies/ml)')
    axs[1][i].set_ylabel("Group Size")
    axs[1][i].get_legend().remove()
    axs[1][i].xaxis.set_tick_params(labelbottom=True)
    fig.align_ylabels([axs[1][i], axs[0][i]])


plt.tight_layout()
plt.savefig(outDir + 'supp_within'+str(NUM_BOOTSTRAPS)+'.jpg', dpi = 300)
plt.close()