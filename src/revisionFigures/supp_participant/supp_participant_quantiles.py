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
PAR_LIST = ['p3', 'p4', 'p7']
EXAMPLE_THRESHOLD = 50000


#For running on desktop
dataDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/"
vlDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/"
outDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/revision/supp_figs/supp_participant/"

# #For running on cluster
# dataDir = "/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/"
# vlDir = "/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/"
# outDir = "/net/feder/vol1/home/evromero/2021_hiv-rec/results/revision/supp_figs/supp_participant/"

#Make the dataframe containing D' ratios
stat_df = zu.combine_drats(dataDir)
stat_df = zu.label_vl_drats(stat_df, vlDir)

stat_df = stat_df[stat_df['Fragment'] != 'F5']
stat_df = stat_df[stat_df['Time_Diff'].gt(50)]

#Only get values for which the viral load doesn't leave the boundaries of the 
#initial and final measurments
stat_df = stat_df[stat_df['Monotonic'] == True]

stat_df = stat_df[stat_df['d_i'] > 0.2]

print(stat_df['Participant'].value_counts())

############################# Panel 1 Analysis ################################
all_par_ests = []
all_group_fits = []
group_size_df = []
x_vals = stat_df['Dist_X_Time'].describe()
par_bins = {}


for curr_par in PAR_LIST:
    #Get the dataframe for only the current participant
    
    curr_stat_df = stat_df[stat_df['Participant'] == curr_par]
    #make a deep copy so we can set values on it
    curr_stat_df = curr_stat_df.copy()


    curr_stat_df['VL_Quantile'], bins = pd.qcut(curr_stat_df['Ave_VL'], q = [0, 0.5, 1], retbins= True, labels = ["Low_VL","High_VL"])

    print('Participant: ' + curr_par)
    print('High D Prime Ratios: ' + str(len(curr_stat_df[curr_stat_df['VL_Quantile'] == 'High_VL'])))
    print('Low D Prime Ratios: ' + str(len(curr_stat_df[curr_stat_df['VL_Quantile'] == 'Low_VL'])))

    #Save the bins for the participant in a dictionary
    par_bins[curr_par] = bins

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
for name, group in all_par_ests.groupby(['Group', 'Participant']):
    lower_conf = np.quantile(group['Estimated_Rho'], 0.025)
    mid_conf = np.quantile(group['Estimated_Rho'], 0.5)
    upper_conf = np.quantile(group['Estimated_Rho'], 0.975)

    #append the confidence interval
    all_conf_ints.append([lower_conf, mid_conf, upper_conf, name[0], name[1]])

all_conf_ints = pd.DataFrame(all_conf_ints, 
    columns=['lower_conf', 'mid_conf', 'upper_conf', 'Group', 'Participant'])


#Output the data that made the figures
all_par_ests.to_csv(outDir + 'all_par_ests_'+ str(NUM_BOOTSTRAPS)+ ' .csv')
all_conf_ints.to_csv(outDir + 'all_conf_ints_' + str(NUM_BOOTSTRAPS)+ '.csv')
print(all_conf_ints)

############################# We want 3 panels ################################
#plot the estimates to show how accurate they are
params = {'figure.figsize':(7, 8), 'axes.labelsize': 6,'axes.titlesize':6,  
          'legend.fontsize': 6, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
          'legend.title_fontsize': 6}
rcParams.update(params)
fig, axs = plt.subplots(len(PAR_LIST),3)


print(all_group_fits)
rcParams.update({'font.size': 8})
linewidth = 1
markersize = 4
markeredgewidth = 1


for i in range(len(PAR_LIST)):
    curr_par = PAR_LIST[i] 
    print(curr_par)
    curr_par_ests = all_par_ests[all_par_ests['Participant'] == curr_par]
    curr_group_fits = all_group_fits[all_group_fits['Participant'] == curr_par]

    ############################ Plot viral loads #################################
    ax0 = axs[i][0]

    #data frame of viral loads for each patient
    viralLoadData = zu.make_viral_load_df(vlDir)
    viralLoadData = viralLoadData[viralLoadData['Participant'] == curr_par]

    #now plot the viral loads
    myplot = sns.lineplot(x = 'Days from infection', y = 'Viral load [virions/ml]',
                           data = viralLoadData, errorbar = None, ax = ax0, color = 'k',
                           markers = True, linewidth = linewidth, markersize = markersize)
    
    color_list = ['tab:blue', 'tab:orange']
    fill_x = np.linspace(0, 6000, 100)
    print(par_bins[curr_par])
    ax0.fill_between(fill_x, 0, par_bins[curr_par][1], color = color_list[0], alpha = 0.3)
    ax0.fill_between(fill_x, par_bins[curr_par][1], 10**6, color = color_list[1], alpha = 0.3)
    ax0.set_ylim(4*10**2, 10**6)
    ax0.set_xlim(0, 6000)
    ax0.set_title('Participant ' + str(curr_par[1:]))
    ax0.set_yscale('log')
    ax0.set_xlabel("Time (days since EDI)")
    ax0.set_ylabel("Viral Load (copies/mL)")


    ########################## Plot the box plots ################################
    ax1 = axs[i][2]

    lower_label = "Below Median \n< " + str(int(par_bins[curr_par][1]))
    upper_label = "Above Median \n> " + str(int(par_bins[curr_par][1]))


    label_dict = {'Low_VL': lower_label, 'High_VL': upper_label}
    curr_par_ests = curr_par_ests.copy()
    curr_par_ests['Group'] = [label_dict[x] for x in curr_par_ests['Group']]

    sns.boxplot(x ='Group', y ='Estimated_Rho', data = curr_par_ests, ax = ax1,  fliersize = 2, linewidth = linewidth,
                palette = 'tab10')

    ax1.set_ylabel(r'Estimated Recombination Rate ($\hat{\rho}$)')
    ax1.set_xlabel(r'Viral Load Quantile (copies/mL)')
    ax1.set_ylim(0, 0.0001)
    ax1.set_title('Participant ' + str(curr_par[1:]))
    ax1.ticklabel_format(style = 'sci', axis = 'y', scilimits = (0,0))
    ax1.xaxis.set_tick_params(labelbottom=True)
   

    ############################# Plotting Panel D ################################
    ax2 = axs[i][1]

    low_50k = curr_group_fits[curr_group_fits['Group'] == 'Low_VL']
    high_50k = curr_group_fits[curr_group_fits['Group'] == 'High_VL']
    sns.lineplot(x ='bin_edges', y ='ratio_bins', data = low_50k, color = 'tab:blue', linewidth = linewidth, label = 'Low VL', ax = ax2)
    sns.lineplot(x ='bin_edges', y ='mid_conf', data = low_50k, color = 'navy',linewidth = linewidth, ax = ax2)


    sns.lineplot(x ='bin_edges', y ='ratio_bins', data = high_50k, color = 'tab:orange', linewidth = linewidth, label = 'High VL', ax = ax2)
    sns.lineplot(x ='bin_edges', y ='mid_conf', data = high_50k, color = 'saddlebrown', linewidth = linewidth, ax = ax2)

    handles, labels = ax2.get_legend_handles_labels()
    labels = ['Below Median', 'Above Median']
    by_label = dict(zip(labels, handles))
    ax2.set_title('Participant ' + str(curr_par[1:]))
    ax2.legend(by_label.values(), by_label.keys(), title = 'Viral Load Group')
    ax2.set_xlabel(r'Distance $\cdot$ Time (bp $\cdot$ generations)')
    ax2.set_ylabel("Linkage Decay Measure (LDM)")

plt.tight_layout()
plt.savefig(outDir + 'supp_participant_' + str(NUM_BOOTSTRAPS) + '.jpg', dpi = 300)
plt.close()