
import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import estimation_util as est_util
import slimUtil as slim_util
from scipy import stats
from matplotlib import rcParams
from matplotlib.lines import Line2D
from scipy.stats import binned_statistic


#This figure is in response to Editor comment 1.1
#First, as you mention in lines 102-106, linkage disequilibrium results from 
#new mutations that occur in particular genetic background and then 
#disassociate with recombination. Then, what you may expect from a within-host
#population of HIV is an equilibrium level of LD due to recurrent mutation and
#recombination, particularly when one anticipates a characteristic genealogy
#(phylogeny) of within-host HIV (for example, Figure 1I of Grenfell et al. 
#2004) that shows continuous turn-over of viral lineages. I therefore expect
#that what is driving the relationship between D' ratio and distance*time in
#Figure 2 is mostly distance (at a given time). I wonder whether you actually
#see D' decrease over time while distance is fixed.

#In this file I am going to plot the linkage curves holding distance constant
#I am going to do this by uniformly sampling distances at each timepoint


NUM_BOOTSTRAPS = 1000
NUM_REPS = 1000
NUM_GROUPS = 10

#The rho values to show in the accuracy panel
ACC_LIST = [r"$2\times10^{-6}$",
                        r"$5\times10^{-6}$",
                        r"$10^{-5}$",
                        r"$2\times10^{-5}$",
                        r"$5\times10^{-5}$",
                        r"$10^{-4}$",
                        r"$2\times10^{-4}$",
                        r"$5\times10^{-4}$",
                        r"$10^{-3}$"]

ACC_LIST = [r"$2\times10^{-5}$", r"$10^{-3}$"]

dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/revision/simulated_estimates_selection/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/revision/reviewer_dist/'



###############################################################################
###############################################################################
###############################################################################
#Load the dataframes that will be plotted
all_stat_dfs = pd.read_pickle(dataDir + "all_stat_dfs_"+ \
                                        str(NUM_BOOTSTRAPS) + "_" + \
                                            str(NUM_GROUPS) + ".pkl")

all_stat_dfs = slim_util.label_rho(all_stat_dfs)
all_stat_dfs['Dist'] = all_stat_dfs['Locus_2'] - all_stat_dfs['Locus_1']


#Get only the D' ratios 50bp apart
all_stat_dfs = all_stat_dfs[all_stat_dfs['Sim_Rho'].isin(ACC_LIST)]
print(all_stat_dfs.head())
print(all_stat_dfs['Time_Diff'].unique())
#Fix the time difference so every point has the correct spacing in the line plot
all_stat_dfs = all_stat_dfs[all_stat_dfs['Time_Diff'] == 150]
all_stat_dfs = all_stat_dfs[all_stat_dfs['Dist'].between(50, 60)]
# all_stat_dfs = all_stat_dfs[all_stat_dfs['Time_1'] == 200]
print(all_stat_dfs.head())
print(len(all_stat_dfs))

#Now plot the linkage trajectories
# num_seg = []
# for curr_rho in ACC_LIST:
#     curr_stat_df = all_stat_dfs[all_stat_dfs['Sim_Rho'] == curr_rho]
#     for name, group in curr_stat_df.groupby(['iter_group', 'rep', 'Time_1', 'Time_2']):
#         num_seg.append(len(set(group['Locus_1'] + group['Locus_2'])))
        
# print(np.mean(num_seg))
# print(np.std(num_seg))
# print(np.median(num_seg))
        


#Now plot the linkage trajectories
for curr_rho in ACC_LIST:
    dist_list = []
    curr_stat_df = all_stat_dfs[all_stat_dfs['Sim_Rho'] == curr_rho]
    my_data = []
    for name, group in curr_stat_df.groupby(['iter_group', 'rep','Locus_1', 'Locus_2']):
        curr_dist = group['Dist'].unique()[0]
        dist_list.append(curr_dist)
        if len(group) == 2:
            group['name'] = name[2] + name[3]
            #Get the first timepoint
            time_1 = min(group['Time_1'])
            d_1 = group[group['Time_1'] == time_1]['d_i'].unique()[0]
            time_2 = max(group['Time_1'])
            d_2 = group[group['Time_1'] == time_2]['d_i'].unique()[0]
            time_3 = max(group['Time_2'])
            d_3 = group[group['Time_1'] == time_2]['d_i'].unique()[0]
            # d_3 = -1 * d_3
            # d_3 = np.exp(-d_3) * d_2
            
            curr_d_vals = [d_1, d_2]
            curr_time_vals = [time_1, time_2]
            order_vals = [1, 2]
            curr_df = pd.DataFrame(list(zip(curr_time_vals, curr_d_vals, order_vals)),columns =['Time', 'd_i', 'order'])
            curr_df['name'] = name[2] + name[3]
            my_data.append(curr_df)
            # sns.lineplot(data=curr_df, x='Time', y='d_i', 
            #              color = 'black', alpha = 0.1, errorbar=None)
    my_data = pd.concat(my_data, ignore_index=True)
    sns.lineplot(data=my_data, x='Time', y='d_i', hue = 'name', errorbar=None)    
    

    timepoint_list = my_data['Time'].unique()
    timepoint_list.sort() 
    print(timepoint_list)

    for i in range(0, len(timepoint_list)-1):
        first_time = timepoint_list[i]
        time_1_data = my_data[my_data['Time'] == first_time]
        time_1_data = time_1_data[time_1_data['order'] == 1]
        d_1 = np.mean(time_1_data['d_i'])

        second_time = timepoint_list[i+1]
        time_2_data = my_data[my_data['Time'] == second_time]
        time_2_data = time_2_data[time_2_data['order'] == 2]
        d_2 = np.mean(time_2_data['d_i'])
        my_aves = pd.DataFrame({'Time': [first_time, second_time], 'd_i': [d_1, d_2]})
        print(my_aves)
        
        sns.lineplot(data=my_aves, x='Time', y='d_i', color = 'red',
                    linewidth = 2, errorbar=None)
    
    redline = Line2D([0], [0], color='r', label = "Average D\' Decay", linewidth= 2)

    plt.legend(handles = [redline], loc = 'upper left', frameon = True)

        




#     for name, group in curr_stat_df.groupby(['iter_group', 'rep']):
#         if len(group) > 1:
#             sns.lineplot(data=group, x='Time_Diff', y='d_ratio', 
#                         color = 'black', alpha = 0.5)
    plt.xlabel(r'Time (generations)')
    plt.ylabel("D\'")
    
    plt.title(r'$\rho$='+curr_rho)
    plt.savefig(outDir + 'dist_fixed_' + curr_rho + '.png', dpi = 300)
    plt.close()

    sns.histplot(dist_list)
    plt.savefig(outDir + 'dist_hist_' + curr_rho + '.png', dpi = 300)
    plt.close()



# #Now plot the linkage trajectories
# for curr_rho in ACC_LIST:
#     curr_stat_df = all_stat_dfs[all_stat_dfs['Sim_Rho'] == curr_rho]
#     for name, group in curr_stat_df.groupby(['Locus_1', 'Locus_2']):
#         print(len(group))
#         if len(group) > 1:
#             print(group)

#     for name, group in curr_stat_df.groupby(['iter_group', 'rep']):
#         if len(group) > 1:
#             sns.lineplot(data=group, x='Time_Diff', y='d_ratio', 
#                         color = 'black', alpha = 0.5)
#     plt.xlabel(r'Time (generations)')
#     plt.ylabel("Linkage Decay Measure")
#     plt.savefig(outDir + 'dist_fixed_' + curr_rho + '.png', dpi = 300)
#     plt.close()

    #Now plot the linkage trajectories of the D' ratios

# ############################# Run the analysis #################################
# all_rho_bins_time = []

# #Bin the points for each rho value
# for curr_rho in all_stat_dfs['Sim_Rho'].unique():
#     if curr_rho not in ACC_LIST:
#         continue
#     curr_stat_df = all_stat_dfs[all_stat_dfs['Sim_Rho'] == curr_rho]

#     #Bin the d' ratios based on time only
#     #Make sure every bin has the same number of points
#     uniform_sampled_dist = []

#     #loop through the possible distances and sample 500 points at each distance
#     for curr_time in curr_stat_df['Time_Diff'].unique():
#         time_lens = []
#         curr_time_df = curr_stat_df[curr_stat_df['Time_Diff'] == curr_time]
#         for i in range(1, 80):
#             curr_dist = curr_time_df[curr_time_df['Dist'] == i]
#             if len(curr_dist) > 30:
#                 uniform_sampled_dist.append(curr_dist.sample(n=30))
#             else:
#                 uniform_sampled_dist.append(curr_dist)
#                 print('***************')
#                 print(i)
#                 print(len(curr_dist))
            
            
#     uniform_sampled_dist = pd.concat(uniform_sampled_dist, ignore_index=True)
    

#     binned_rat_time, binedges_time, bin_nums_time = binned_statistic(
#         np.array(uniform_sampled_dist['Time_Diff'].tolist()), 
#         np.array(uniform_sampled_dist['d_ratio'].tolist()), bins = 20)
    
    
#     curr_rho_bins_time = pd.DataFrame({'Time': binedges_time[1:], 
#                                        'D_Ratio': binned_rat_time})
#     curr_rho_bins_time['Sim_Rho'] = curr_rho
#     curr_rho_bins_time['Sim_float_rho'] = uniform_sampled_dist['Sim_float_rho'].unique()[0]

#     all_rho_bins_time.append(curr_rho_bins_time)


# all_rho_bins_time = pd.concat(all_rho_bins_time, ignore_index=True)


# all_rho_bins_time['Sim_Rho'] = all_rho_bins_time['Sim_Rho'].astype('category')

# ######################### Configure Plot Settings #############################
# #plot the estimates to show how accurate they are
# params = {'figure.figsize':(7.5, 3), 'axes.labelsize': 6,'axes.titlesize':6,  
#           'legend.fontsize': 6, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
#           'legend.title_fontsize': 6}
# linewidth = 1
# markerSize = 3
# rcParams.update(params)
# fig, axs = plt.subplots(1, 2)

# ########################## Panel C Time #######################################
# sns.lineplot(data=all_rho_bins_time, x='Time', y='D_Ratio', ax = axs[1],
#              hue='Sim_Rho', linewidth = 1, 
#              palette=sns.color_palette("viridis", 
#              n_colors=len(all_rho_bins_time['Sim_Rho'].unique())), 
#              hue_order = ACC_LIST)

# axs[1].set_xlabel(r'Time (generations)')
# axs[1].set_ylabel("Linkage Decay Measure")
# # axs[1].get_legend().remove()


# sns.histplot(uniform_sampled_dist, x='Dist', hue='Time_Diff', alpha = 0.3,
#              ax = axs[0], palette=sns.color_palette("YlOrBr", 
#              n_colors=len(uniform_sampled_dist['Time_Diff'].unique())))
# plt.savefig(outDir + 'dist_fixed.png', dpi = 300)
# plt.close()


# # plt.savefig(outDir + 'time_fixed_150.png', dpi = 300)
