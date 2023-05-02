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
outDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/paper/supp_loo/"

#For running on cluster
# dataDir = "/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/"
# vlDir = "/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/"
# outDir = "/net/feder/vol1/home/evromero/2021_hiv-rec/results/paper/supp_loo/"

#Make the dataframe containing D' ratios
stat_df = zu.combine_drats(dataDir)
stat_df = zu.label_vl_drats(stat_df, vlDir)

stat_df = stat_df[stat_df['Fragment'] != 'F5']
stat_df = stat_df[stat_df['Time_Diff'].gt(50)]

#Only get values for which the viral load doesn't leave the boundaries of the 
#initial and final measurments
stat_df = stat_df[stat_df['Monotonic'] == True]

all_par_ests = []
all_group_fits = []
group_size_df = []

#Estimate rates specifically excluding each individual
for curr_par in stat_df['Participant'].unique():
    #Get the dataframe for everyone except the current participant
    excl_par_df = stat_df[stat_df['Participant'] != curr_par]
    #Make a deep copy so we can assign some values (the viral load groups)
    excl_par_df = excl_par_df.copy()

    lower_label = "Lower \n<" + str(int(np.quantile(excl_par_df['Ave_VL'], 0.33)))
    middle_label = "Middle \n" + str(int(np.quantile(excl_par_df['Ave_VL'], 0.33))) + " - " + str(int(np.quantile(excl_par_df['Ave_VL'], 0.66))) 
    upper_label = "Upper \n>" + str(int(np.quantile(excl_par_df['Ave_VL'], 0.66)))

    #First we need to group the data by the viral load quantiles
    excl_par_df['VL_Quantile'], bins = pd.qcut(excl_par_df['Ave_VL'], q = QUANTILE_LIST, retbins= True, labels = [lower_label, middle_label, upper_label])
    grouped_stats = excl_par_df.groupby(by = ['VL_Quantile'])

    for name, group in grouped_stats:
        #Get the current estimate
        lower_fit, upper_fit, estimate_df = plne.bootstrap_rho(group,
                                                                NUM_BOOTSTRAPS)
        estimate_df['Group'] = name
        estimate_df['excl_par'] = curr_par

        #Add the size of the group to the dictionary which labels the axis
        group_size_df.append([name, group.shape[0], curr_par])
        
        #Add all of the current results
        all_par_ests.append(estimate_df)


all_par_ests = pd.concat(all_par_ests, ignore_index=True)
all_group_sizes = pd.DataFrame(group_size_df, columns=['Group', 'Size', 'Excl Par'])
all_conf_ints = []

#Calculate the confidence intervals for the estimates with low + high VL
for name, group in all_par_ests.groupby(['Group', 'excl_par']):
    lower_conf = np.quantile(group['Estimated_Rho'], 0.025)
    mid_conf = np.quantile(group['Estimated_Rho'], 0.5)
    upper_conf = np.quantile(group['Estimated_Rho'], 0.975)

    #append the confidence interval
    all_conf_ints.append([lower_conf, mid_conf, upper_conf, name[0], name[1]])

all_conf_ints = pd.DataFrame(all_conf_ints, 
    columns=['lower_conf', 'mid_conf', 'upper_conf', 'Group', 'excl_par'])

#Output the data that made the figures
all_par_ests.to_csv(outDir + 'all_par_ests_quantiles_' + str(NUM_BOOTSTRAPS) + '.csv')
all_conf_ints.to_csv(outDir + 'all_conf_ints_quantiles' + str(NUM_BOOTSTRAPS) + '.csv')

######################### Configure Plot Settings #############################
#plot the estimates to show how accurate they are
params = {'figure.figsize':(7, 8), 'axes.labelsize': 6,'axes.titlesize':6,  
          'legend.fontsize': 6, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
          'legend.title_fontsize': 6}
rcParams.update(params)

#Make a three panel figure
linewidth = 1
markersize = 4

###############################################################################

fig, axs = plt.subplots(4, 3, sharey = True)

ax_list = axs.flatten()
par_list = ['p1', 'p2', 'p3','p4', 'p5', 'p6', 'p7', 'p8', 'p9', 'p10', 'p11']

#plot each of the participant subplots separately
for i in range(len(par_list)):
    curr_excl = par_list[i]
    curr_df = all_par_ests[all_par_ests['excl_par'] == curr_excl]
    curr_ax = ax_list[i]
    sns.boxplot(x ='Group', y ='Estimated_Rho',data = curr_df, ax = curr_ax,
                     fliersize = 2, linewidth = linewidth,
                    palette=['tab:blue', 'tab:gray', 'tab:orange'])
    curr_ax.set_ylabel(r'')
    curr_ax.set_xlabel(r'')
    curr_ax.set_title('Excluding Participant ' + str(curr_excl[1:]), pad = -14)
    curr_ax.ticklabel_format(style = 'sci', axis = 'y', scilimits = (0,0))
    curr_ax.xaxis.set_tick_params(labelbottom=True)

axs[0][0].set_ylabel(r'Estimated Recombination Rate ($\hat{\rho}$)')
axs[1][0].set_ylabel(r'Estimated Recombination Rate ($\hat{\rho}$)')
axs[2][0].set_ylabel(r'Estimated Recombination Rate ($\hat{\rho}$)')
axs[3][0].set_ylabel(r'Estimated Recombination Rate ($\hat{\rho}$)')
axs[3][0].set_xlabel(r'Viral Load Tertile (copies/mL)')
axs[3][1].set_xlabel(r'Viral Load Tertile (copies/mL)')
axs[2][2].set_xlabel(r'Viral Load Tertile (copies/mL)')

plt.subplots_adjust(hspace= 0.4)
plt.savefig(outDir + "loo_estimates_quant_" + str(NUM_BOOTSTRAPS) + ".jpg", dpi = 300)
plt.close()