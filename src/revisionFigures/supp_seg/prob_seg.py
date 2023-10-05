import sys
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import zaniniUtil as zu
import slimUtil as slim_util
from matplotlib import rcParams
from matplotlib.lines import Line2D
from scipy.stats import binned_statistic

NUM_BOOTSTRAPS = 1000
NUM_GROUPS = 10

#In this file I am going to match the amount of data in each viral load tertile as closely as possible
# to the amount of data in each selection simulation group.
simDataDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/revision/simulated_estimates_selection/"

invivoDataDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/"
vlDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/"
outDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/revision/supp_figs/supp_seg/"

params = {'figure.figsize':(4, 2.5), 'axes.labelsize': 6,'axes.titlesize':6,  
          'legend.fontsize': 6, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
          'legend.title_fontsize': 6}
rcParams.update(params)

###############################################################################
############################## Simulated Data #################################
###############################################################################
all_stat_dfs_sim = pd.read_pickle(simDataDir + "all_stat_dfs_"+ \
                                        str(NUM_BOOTSTRAPS) + "_" + \
                                            str(NUM_GROUPS) + ".pkl")
all_stat_dfs_sim = slim_util.label_rho(all_stat_dfs_sim)


sim_loci_counts = []
#Count the number of segregating sites per rep
for name, group in all_stat_dfs_sim.groupby(['rep', 'Sim_Rho', 'Time_1']):
    curr_loci = set(group['Locus_1'].unique().tolist() + \
                    group['Locus_2'].unique().tolist())
    iter_group = group['iter_group'].unique().tolist()[0]
    num_drat = len(group)
    sim_loci_counts.append([name[0], name[1], name[2],# name[3],
                            iter_group, len(curr_loci), num_drat])

sim_loci_counts = pd.DataFrame(sim_loci_counts,
                    columns = ['rep', 'Sim_Rho', 'Time_1',
                             'iter_group', 'Segregating Sites', 'LDMs'])
sim_loci_counts['Segregating Sites'] = sim_loci_counts['Segregating Sites']/700

print('The total number of segregating sites is')
print(sim_loci_counts.groupby(['Sim_Rho']).sum()['Segregating Sites']*700)
print("The median number of segregating sites per timepoint pair is")
print(sim_loci_counts.groupby(['Sim_Rho']).median()['Segregating Sites'])
print("The mean number of segregating sites per timepoint pair is")
print(sim_loci_counts.groupby(['Sim_Rho']).mean()['Segregating Sites'])

###############################################################################
################################ In Vivo Data #################################
###############################################################################
fragmentLen = {
    'F1': 1823,
    'F2': 1894,
    'F3': 1583,
    'F4': 1696,
    'F6': 1829
}

#Make the dataframe containing D' ratios
iv_df = zu.combine_drats(invivoDataDir)
iv_df = iv_df[iv_df['Fragment'] != 'F5']
iv_df = iv_df[iv_df['Time_Diff'].gt(50)]
# iv_df = iv_df[iv_df['Participant'] != 'p4']
# iv_df = iv_df[iv_df['Participant'] != 'p7']
iv_df = iv_df[iv_df['d_i'] > 0.2]
iv_df = iv_df[iv_df['Dist_X_Time'] < 50000]
iv_df = zu.label_vl_drats(iv_df, vlDir)

#Summarize the D' ratios and segregating loci per participant
iv_df = iv_df.groupby(['Participant', 'Fragment', 'Time_1'])
summary_df = []
for name, curr_group in iv_df:
    viral_load = (curr_group['Ave_VL'].unique())[0]
    num_seg = curr_group['Locus_1'] + curr_group['Locus_2']
    num_seg = len(num_seg.unique())
    num_drat = len(curr_group)
    summary_df.append([name[0], name[1], name[2],# name[3],
                       num_seg, num_drat, viral_load])

inv_loci_counts = pd.DataFrame(summary_df, columns=['Participant', 'Fragment',
                'Time_1', 'Segregating Sites',  'LDMs',
                'Viral Load'])
inv_loci_counts['Viral Load'] = np.log10(inv_loci_counts['Viral Load'])
frag_lens = []

for index,row in inv_loci_counts.iterrows():
    curr_frag = row['Fragment']
    frag_lens.append(fragmentLen[curr_frag])
inv_loci_counts['Fragment Length'] = frag_lens
inv_loci_counts['Segregating Sites'] = inv_loci_counts['Segregating Sites'] / \
    inv_loci_counts['Fragment Length']

###############################################################################
############################## Plot the Results ###############################
###############################################################################

sns.scatterplot(data= inv_loci_counts, x='Time_1', label = 'One In Vivo Sequencing Fragment',
                y='Segregating Sites', color = 'k', alpha = 0.3, s = 7)
sns.scatterplot(data= sim_loci_counts, x='Time_1', label = 'One Simulation',
                y='Segregating Sites', color = 'red', alpha = 0.3, s = 7)

#plot binned averages of the results
my_bins = range(0, 1200, 200)

bin_count, binedges, bin_nums = binned_statistic(
        sim_loci_counts['Time_1'].to_numpy(), 
        sim_loci_counts['Segregating Sites'].to_numpy(),
        bins = my_bins, statistic= 'median')

sns.lineplot(x = binedges[:-1], y = bin_count, color = 'red', label = 'Simulated Median')

bin_count, binedges, bin_nums = binned_statistic(
        inv_loci_counts['Time_1'].to_numpy(), 
        inv_loci_counts['Segregating Sites'].to_numpy(),
        bins = my_bins, statistic= 'median')

sns.lineplot(x = binedges[1:], y = bin_count, color = 'k', label = 'In Vivo Median')
plt.ylabel('Proportion of Segregating Sites')
plt.xlabel('Time Point (generations)')
plt.tight_layout()
plt.savefig(outDir + 'segregating_sites_median_single_time.png', dpi = 300)
plt.close()