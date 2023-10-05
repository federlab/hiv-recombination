import os
import sys
import pandas as pd
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import seaborn as sns
import numpy as np
import estimation_util as est_util
import slimUtil as slim_util
import zaniniUtil as zu
import plot_neher as plne
import matplotlib.pyplot as plt
import autocorrelation as autocorr
from scipy import stats
from matplotlib.lines import Line2D
from matplotlib import rcParams

#In this file I am going to be making a figure that illustrates the effects
#of decreased sampling depth
dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/revision/simulated_estimates_sampling/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/revision/supp_figs/supp_sampling/'
inv_dataDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/"
vlDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/"

#For running on cluster
# inv_dataDir = "/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/"
# vlDir = "/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/"
# dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/revision/simulated_estimates_sampling/'
# outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/supp_figs/supp_sampling/'

###############################################################################
##################### Setting the figure parameters ###########################
###############################################################################
#plot the estimates to show how accurate they are
params = {'figure.figsize':(7, 3), 'axes.labelsize': 6,'axes.titlesize':6,  
          'legend.fontsize': 6, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
          'legend.title_fontsize': 6}
rcParams.update(params)

fig, axs = plt.subplots(1, 2, sharey= True, gridspec_kw={'width_ratios':[2.5,1]})
plt.subplots_adjust(wspace = 0.05, hspace=0)
ax0 = axs[0]
ax1 = axs[1]

linewidth = 1

###############################################################################
######################## Plotting the simulated estimates #####################
###############################################################################

DIST_TIME_MAX = 50000
NUM_BOOTSTRAPS = 1000
NUM_REPS = 5000
NUM_GROUPS = 5

SAMPLE_RHO_LIST = [ r"$10^{-4}$", r"$10^{-5}$"]
SAMPLING_DEPTH = ['10', '20', '40', '80']

order_list = []
for depth in SAMPLING_DEPTH:
    for rho in SAMPLE_RHO_LIST:
        my_label = r'$\rho$ = ' + rho + ' \n Depth = ' + depth
        order_list.append(my_label)

all_stat_dfs = pd.read_pickle(dataDir + "all_stat_dfs_"+ \
                                        str(NUM_BOOTSTRAPS) + "_" + \
                                            str(NUM_GROUPS) + ".pkl")

# print(all_stat_dfs.head())
for depth in all_stat_dfs['depth'].unique():
    print('***************************************')
    print('Depth is '+ str(depth))
    print(slim_util.count_seg_sites(all_stat_dfs[all_stat_dfs['depth'] == depth]))


all_ests_df = pd.read_pickle(dataDir+ "all_estimate_df_"+ \
                                        str(NUM_BOOTSTRAPS) + "_" + \
                                            str(NUM_GROUPS) + ".pkl")


#Relabel the Rho values
all_ests_df = slim_util.label_rho(all_ests_df)
all_ests_df = all_ests_df[all_ests_df['Sim_Rho'].isin(SAMPLE_RHO_LIST)]

#Make a column that gives the rho and sampling info
all_ests_df['Rho_Sampling'] = r'$\rho$ = ' + all_ests_df['Sim_Rho'] +\
                     ' \n Depth = ' + all_ests_df['depth'].astype(str)

my_palette = sns.color_palette(['k'])
sns.stripplot(x = 'Rho_Sampling', y = 'Estimated_Rho', hue = 'Group',
                data = all_ests_df, jitter = True, s = 2, palette = my_palette,
                alpha = 0.05, order = order_list, dodge = True, ax = ax0)

# distance across the "X" or "Y" stipplot column to span, in this case 40%
label_width = 0.5
linewidth = 1

# we need to save the figure once here so that we can access the axis labels
plt.savefig(outDir + "sampling_accuracy_bootstraps_strip.jpg", dpi = 300)
            
for tick, text in zip(ax0.get_xticks(), ax0.get_xticklabels()):
    sample_name = text.get_text()  # "X" or "Y"

    #get the float value of rho corresponding with the tick
    rho_val = all_ests_df[all_ests_df['Rho_Sampling'] == sample_name]
    rho_val = rho_val['Sim_float_rho'].unique()[0]
    estimate = np.mean(all_ests_df[all_ests_df['Rho_Sampling'] == sample_name]['Estimated_Rho'])

    #plot horizontal lines across the column, centered on the tick
    ax0.plot([tick-label_width/2, tick+label_width/2], [rho_val, rho_val],
            lw=1, color='r', linestyle = '--')

#make a legend showing which line is the mean and which line is the truth
redline = Line2D([0], [0], color='r', linestyle = '--', label = r'True $\rho$', linewidth= linewidth)

ax0.set_ylim(0.00000001, 0.01)
ax0.set_xticklabels(ax0.get_xticklabels(), rotation = 90)
ax0.set_xlabel(r'True Recombination Rate ($\rho$) and Sampling Depth')
ax0.set_ylabel(r'Estimated Recombination Rate $(\hat{\rho})$')
ax0.set_yscale('log')
ax0.get_legend().remove()
ax0.legend(handles = [redline], loc = 'upper left', frameon = True)


###############################################################################
######################### Plotting the Zanini estimates #######################
###############################################################################
NUM_BOOTSTRAPS = 1000
EXAMPLE_THRESHOLD = 50000
QUANTILE_LIST = [0, 0.33, 0.66, 1]


#Make the dataframe containing D' ratios
stat_df = zu.combine_drats(inv_dataDir)
stat_df = zu.label_vl_drats(stat_df, vlDir)

stat_df = stat_df[stat_df['Fragment'] != 'F5']
stat_df = stat_df[stat_df['Time_Diff'].gt(50)]

#Only get values for which the viral load doesn't leave the boundaries of the 
#initial and final measurments
stat_df = stat_df[stat_df['Monotonic'] == True]

all_par_ests = []
all_group_fits = []
group_size_df = []

lower_label = "Zanini et al.\n Lower Tertile\n<" + str(int(np.quantile(stat_df['Ave_VL'], 0.33)))
middle_label = "Zanini et al.\n Middle Tertile\n" + str(int(np.quantile(stat_df['Ave_VL'], 0.33))) \
                + " - " + str(int(np.quantile(stat_df['Ave_VL'], 0.66))) 
upper_label = "Zanini et al.\n Upper Tertile\n>" + str(int(np.quantile(stat_df['Ave_VL'], 0.66)))

#First we need to group the data by the viral load quantiles
stat_df['VL_Quantile'], bins = pd.qcut(stat_df['Ave_VL'], q = QUANTILE_LIST, 
            retbins= True, labels = [lower_label, middle_label, upper_label])
grouped_stats = stat_df.groupby(by = ['VL_Quantile'])


for name, group in grouped_stats:
    #Get the current estimate
    lower_fit, upper_fit, estimate_df = plne.bootstrap_rho(group,
                                                            NUM_BOOTSTRAPS)
    estimate_df['Group'] = name

    #Add the size of the group to the dictionary which labels the axis
    group_size_df.append([name, group.shape[0]])
    
    #Add all of the current results
    all_par_ests.append(estimate_df)
    # all_group_fits.append(group_fits)

all_par_ests = pd.concat(all_par_ests, ignore_index=True)
# all_group_fits = pd.concat(all_group_fits, ignore_index=True)
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

print(all_par_ests.head())

###############################################################################
######################### Plotting the Caskey estimates #######################
###############################################################################
caskey_data = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/caskey/'

all_stat_df_caskey = []
total_loci = 0

#Loop through each of the directory files and calculate the D' ratios
for curr_dir in os.listdir(caskey_data):
    if curr_dir[0] == '.':
        continue
    if not os.path.isdir(caskey_data + curr_dir):
        continue

    inFile = caskey_data + curr_dir + '/linkage/r2_and_D'

    #Get the D' ratios
    stat_df = autocorr.calculate_d_ratios(inFile)
    # stat_df = stat_df[stat_df['Time_2'] <= 14]
    stat_df = stat_df[stat_df['Dist_X_Time'] <= 50000]
    if len(stat_df) == 0:
        continue

    all_stat_df_caskey.append(stat_df)

all_stat_df_caskey = pd.concat(all_stat_df_caskey, ignore_index=True)

lower_fit, upper_fit, estimate_df = plne.bootstrap_rho(all_stat_df_caskey,
                                                            NUM_BOOTSTRAPS)
estimate_df['Group'] = 'Caskey et al. \n mean = 12766\n sd = 9859'

print('The median Caskey estimate is')
print(np.median(estimate_df['Estimated_Rho']))
print('The 95% confidence interval is')
print(np.quantile(estimate_df['Estimated_Rho'], 0.025))
print(np.quantile(estimate_df['Estimated_Rho'], 0.975))

all_par_ests = pd.concat([estimate_df, all_par_ests], ignore_index=True)

#Now estimate 


my_palette=['tab:red', 'tab:blue', 'tab:gray', 'tab:orange']
sns.stripplot(x ='Group', y = 'Estimated_Rho', hue = 'Group',
                data = all_par_ests, jitter = 1/40, s = 2, palette = my_palette,
                alpha = 0.05, ax = ax1)
            
ax1.set_xticklabels(ax1.get_xticklabels(), rotation = 90)
ax1.set_xlabel('Data Set and Viral Load (copies/mL)')
ax1.set_ylabel(r'Estimated Recombination Rate $(\hat{\rho})$')
ax1.set_yscale('log')
ax1.get_legend().remove()

plt.tight_layout()

print(outDir + "sampling_accuracy_bootstraps_strip.jpg")
plt.savefig(outDir + "sampling_accuracy_bootstraps_strip.jpg", dpi = 300)