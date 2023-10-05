import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import slimUtil as slim_util
from matplotlib import rcParams

#The length of the simulated segments
GENOME_LENGTH = 700
NE_LIST = ['5000', '10000', '50000']

#Today I am looking through varying NE simulation data to see if the number of 
#segregating sites, sequence divergence etc match the in vivo data

sim_dataDir = "/Volumes/feder-vol1/project/hiv_recombination/data/simulated/2023_08_28_selection_vNE/"
inv_dataFile = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini/elife-11282-fig5-data1-v2/divdiv_AB.tsv" 
outDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/revision/supp_figs/supp_match/"

###############################################################################
############################### Plot Parameters ###############################
params = {'figure.figsize': (7.5, 4), 'axes.labelsize': 6,'axes.titlesize':6,  
          'legend.fontsize': 6, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
          'legend.title_fontsize': 6}
linewidth = 1
markerSize = 3
rcParams.update(params)


###############################################################################
############################ Simulated Results ################################
###############################################################################
all_dist_df = []

#Loop through the simulations
for curr_dir in os.listdir(sim_dataDir):
    if not os.path.isdir(sim_dataDir+curr_dir) or curr_dir[0] == '.':
        continue
    
    
    #Collect the info for the run
    run_info = curr_dir.split('_')


    rho_val = run_info[1]
    rho_val = rho_val[3:]
    ne_val = run_info[2]
    ne_val = ne_val[2:]
    if ne_val not in NE_LIST:
        continue
    rep = run_info[-1]
    rep = rep[3:]

    #For now only look at the first 10 reps
    if int(rep) > 1:
        continue


    curr_file = sim_dataDir + curr_dir + "/slim_output.txt"
    dist_df = slim_util.calc_pop_metrics(curr_file, GENOME_LENGTH)

    #Convert to years
    dist_df['years'] = (dist_df['timepoint'] * 2)/365
    dist_df['Ne'] = ne_val
    dist_df['Sim_float_rho'] = float(rho_val)
    dist_df['Sim_Rho'] = rho_val
    all_dist_df.append(dist_df)


all_dist_df = pd.concat(all_dist_df, ignore_index=True)


ave_dist = all_dist_df.groupby(['years', 'Ne', 'Sim_Rho', 'Sim_float_rho']).mean()
ave_dist = ave_dist.reset_index()
ave_dist = slim_util.label_rho(ave_dist)
ave_dist = slim_util.label_Ne(ave_dist)


#Lastly, plot the results
fig, ax = plt.subplots(1, 2, sharex=True)
plt.subplots_adjust(wspace = 0.5)
sns.lineplot(data = ave_dist, x = 'years', y = 'mean_pairwise_dist',
              palette = sns.color_palette("YlGnBu", n_colors=len(NE_LIST)),
              hue = 'Ne', style = 'Sim_Rho', ax = ax[1], 
              style_order = [r"$10^{-3}$", r"$10^{-4}$", 
                             r"$2\times10^{-4}$", r"$10^{-5}$",
                             r"$2\times10^{-5}$", r"$2\times10^{-6}$"],
             hue_order = [r"$N_e=5\times10^3$", r"$N_e=10^4$", r"$N_e=5\times10^4$"])

ax[1].set_xlabel('Years')
ax[1].set_ylim(0, 0.03)
ax[1].set_ylabel('Mean Pairwise Hamming Distance')
                    
#Now I will measure sequence divergence (average hamming distance from founder)
sns.lineplot(data = ave_dist, x = 'years', y = 'mean_sequence_div',
             palette = sns.color_palette("YlGnBu", n_colors=len(NE_LIST)),
             hue = 'Ne', style = 'Sim_Rho', ax = ax[0],
             style_order = [r"$10^{-3}$", r"$10^{-4}$", 
                             r"$2\times10^{-4}$", r"$10^{-5}$",
                             r"$2\times10^{-5}$", r"$2\times10^{-6}$"],
             hue_order = [r"$N_e=5\times10^3$", r"$N_e=10^4$", r"$N_e=5\times10^4$"])
ax[0].set_xlabel('Years')
ax[0].set_ylim(0, 0.045)
ax[0].set_ylabel('Mean Sequence Divergence')


###############################################################################
############################## In Vivo Results ################################
###############################################################################
#First, read in the data from Zanini et al. 2015 Figure 5
#https://doi.org/10.7554/eLife.11282 

#The following lines just reconstruct the timebins used for Zanini et al. 2015 Figure 5
#I am using the code from
#https://github.com/neherlab/HIVEVO_figures/blob/master/src/syn_nonsyn_divdiv.py
#Starting at line 178

time_bins = np.array([0, 200, 500, 1000, 1500, 2000, 4000])
time_binc = 0.5*(time_bins[:-1]+time_bins[1:])
time_binc = time_binc/365.25

#Now read in the data and put the reconstructed labels on it
inv_data = pd.read_csv(inv_dataFile, sep='\t', header=None, index_col=False, 
                        names = ['Metric', 'Type', 'Protein', 
                                 'Time_0', 'Time_1', 'Time_2',
                                 'Time_3', 'Time_4', 'Time_5'])
inv_data = inv_data.melt(id_vars= ['Metric', 'Type', 'Protein'], var_name='Timepoint')

inv_data['Time_x'] = inv_data['Timepoint'].str.replace('Time_', '')
inv_data['Time_x'] = inv_data['Time_x'].astype(int)
inv_data['Time_x'] = time_binc[inv_data['Time_x']]

#Summarize the data based on timepoint and metric
inv_data = inv_data.groupby(['Metric', 'Time_x']).mean()
inv_data = inv_data.reset_index()

inv_diversity = inv_data[inv_data['Metric'] == 'diversity']
inv_divergence = inv_data[inv_data['Metric'] == 'divergence']

#Now plot the results
sns.lineplot(data = inv_diversity, x = 'Time_x', y = 'value', color = 'red', label = 
             'In Vivo', ax = ax[1])

sns.lineplot(data = inv_diversity, x = 'Time_x', y = 'value', color = 'red', ax = ax[0])

print(ax[1].get_legend_handles_labels()[1])

new_labels = ax[1].get_legend_handles_labels()[1]
new_labels[0] = r'$N_e$'
new_labels[4] = r'Simulated $\rho$'

ax[1].legend(loc='upper center', bbox_to_anchor=(-0.45, 1.15), ncol=6,
           labels = new_labels, handles = ax[1].get_legend_handles_labels()[0], frameon = True)

ax[0].get_legend().remove()
# ax[1].get_legend().remove()



plt.xlim(0,8)
# plt.tight_layout
plt.savefig(outDir + 'div_div.png', dpi = 300)