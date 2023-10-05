import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import slimUtil as slim_util
from scipy.stats import binned_statistic
import plot_neher as plne
from scipy import optimize
from sklearn.metrics import mean_squared_error
from matplotlib import rcParams

THRESHOLD = 0.2
DIST_EXAMPLE = 50000
NUM_BOOTSTRAPS = 1000
NUM_STORED_GROUPS = 40

D_PRIME_NUMS = [2500, 5000, 10000, 25000, 50000]
SEG_LOCI_NUMS = [100, 200, 500, 1000, 2000]
DIST_TIME_MAX = [10000, 25000, 50000, 100000, 200000]

NUM_REPS = 200
NUM_RMSE_GROUPS = 10
#Used for the distance cutoff plot
RATIOS_PER_GROUP = 25000

outTag = "_revised.png"

#For running on Cluster
dataDir = "/net/feder/vol1/home/evromero/2021_hiv-rec/results/paper/simulated_estimates/"
outDir = "/net/feder/vol1/home/evromero/2021_hiv-rec/results/paper/supp_params_revised/"

# dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/paper/simulated_estimates/'
# outDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/paper/supp_params_revised/"

#First, we will get all of the data and divide it into groups
all_stat_dfs = pd.read_pickle(dataDir + "all_stat_dfs_"+ \
                                        str(NUM_BOOTSTRAPS) + "_" + \
                                            str(NUM_STORED_GROUPS) + ".pkl")


#Filter the groups so the first d' value is greater than 0.2
all_stat_dfs = all_stat_dfs[all_stat_dfs['d_i'] > THRESHOLD]
all_stat_dfs = slim_util.label_rho(all_stat_dfs)


#We need to redo the grouping because we only want 10 groups that are doublesized
#Randomly divide the reps into 10 groups
rep_groups = np.array(range(0, NUM_REPS+1))
np.random.shuffle(rep_groups)
rep_groups = np.array_split(rep_groups, NUM_RMSE_GROUPS)


#Make a dictionary to label each group
group_dict = {}
for i in range(len(rep_groups)):
    for j in rep_groups[i]:
        group_dict[j] = i

#Add the group labels to the dataframe      
group_labels = [group_dict[x] for x in all_stat_dfs['rep']]
all_stat_dfs['iter_group'] = group_labels

#loop through each of the sample sizes for segregating loci
estimate_df_seg = []
for curr_size in SEG_LOCI_NUMS:
    print("# Seg Loci is: " + str(curr_size))
    #loop through each rho value
    for curr_rho in all_stat_dfs['Sim_Rho'].unique():
        print(curr_rho)
        curr_float_rho = all_stat_dfs[all_stat_dfs['Sim_Rho'] == curr_rho]['Sim_float_rho'].unique()[0]
        curr_rho_stat = all_stat_dfs[all_stat_dfs['Sim_Rho'] == curr_rho]

        #Add the distance example cutoff
        curr_rho_stat = curr_rho_stat[curr_rho_stat['Dist_X_Time'] <= DIST_EXAMPLE] 

        #loop through each iter group
        for curr_iteration in range(1, NUM_RMSE_GROUPS):
            #get the data for the current rho and iteration
            curr_stat_df = curr_rho_stat[curr_rho_stat['iter_group'] == curr_iteration]
            reps_in_group = len(curr_stat_df['rep'].unique())

            sample_list = []
            
            #Next we need to sample segregating sites but we have to group by locus and simulation
            #To do this we'll make columns that give the rep and the locus
            #basically sample an even number of segregating loci from each rep
            for curr_rep in curr_stat_df['rep'].unique():
                rep_stat_df = curr_stat_df[curr_stat_df['rep'] == curr_rep]
                curr_loci = set(rep_stat_df['Locus_1'].unique().tolist() + rep_stat_df['Locus_2'].unique().tolist())

                #now sample a number of segregating loci from the current rep
                curr_sample = np.random.choice(list(curr_loci), size = int(curr_size/reps_in_group), replace = False)
                sample_list.append(rep_stat_df[rep_stat_df['Locus_1'].isin(curr_sample) & rep_stat_df['Locus_2'].isin(curr_sample)])

            
            curr_stat_df = pd.concat(sample_list, ignore_index=True)

            #get the estimate and fit for the current dataset and sample size
            x_vals = curr_stat_df['Dist_X_Time'].unique()
            #Get the current estimate
            lower_fit, upper_fit, estimate_df = plne.bootstrap_rho(curr_stat_df,
                                                                NUM_BOOTSTRAPS)

            lower_conf = np.quantile(estimate_df['Estimated_Rho'], 0.025)
            mid_est = np.quantile(estimate_df['Estimated_Rho'], 0.5)
            upper_conf = np.quantile(estimate_df['Estimated_Rho'], 0.975)

            estimate_df_seg.append([
                mid_est, curr_rho, curr_iteration, 
                curr_size, curr_float_rho])  

estimate_df_seg = pd.DataFrame(estimate_df_seg, columns=["Est_Rho", 
                    'Sim_Rho', 'iter_name', 'sample_size', 'Sim_float_rho'] )


#loop through each of the sample sizes
estimate_df_size = []
for curr_size in D_PRIME_NUMS:
    print("Size is: " + str(curr_size))
    #loop through each rho value
    for curr_rho in all_stat_dfs['Sim_Rho'].unique():
        print(curr_rho)
        curr_rho_stat = all_stat_dfs[all_stat_dfs['Sim_Rho'] == curr_rho]
        curr_float_rho = all_stat_dfs[all_stat_dfs['Sim_Rho'] == curr_rho]['Sim_float_rho'].unique()[0]
        
        #Add the distance example cutoff
        curr_rho_stat = curr_rho_stat[curr_rho_stat['Dist_X_Time'] <= DIST_EXAMPLE] 

        #loop through each iter group
        for curr_iteration in range(1, NUM_RMSE_GROUPS):
            #get the data for the current rho and iteration
            curr_stat_df = curr_rho_stat[curr_rho_stat['iter_group'] == curr_iteration]
            
            #sample segregating sites while the D' num is below the threshold
            sample_df = []
            sample_size = 0
            sampled_loci = set()
            unsampled_loci = set(list(curr_stat_df['Locus_1'].unique()) + list(curr_stat_df['Locus_2'].unique().tolist()))
            while sample_size < curr_size and len(unsampled_loci) > 0:
                #sample another locus
                my_sample = np.random.choice(list(unsampled_loci))
                sampled_loci.add(my_sample)
                unsampled_loci.remove(my_sample)
                sample_df = curr_stat_df[curr_stat_df['Locus_1'].isin(sampled_loci) & (curr_stat_df['Locus_2'].isin(sampled_loci))]
                sample_size = len(sample_df)
            
            # print(sample_size)
            curr_stat_df= sample_df
            
            if len(unsampled_loci) == 0:
                print("No more loci to sample")
                print(curr_iteration)
                quit()

            #get the estimate and fit for the current dataset and sample size
            x_vals = curr_stat_df['Dist_X_Time'].unique()
            #Get the current estimate
            lower_fit, upper_fit, estimate_df = plne.bootstrap_rho(curr_stat_df,
                                                                NUM_BOOTSTRAPS)

            lower_conf = np.quantile(estimate_df['Estimated_Rho'], 0.025)
            mid_est = np.quantile(estimate_df['Estimated_Rho'], 0.5)
            upper_conf = np.quantile(estimate_df['Estimated_Rho'], 0.975)

            estimate_df_size.append([
                mid_est, curr_rho, curr_iteration, curr_size,
                curr_float_rho])  
            

estimate_df_size = pd.DataFrame(estimate_df_size, columns=["Est_Rho", 
                    'Sim_Rho', 'iter_name', 'sample_size', 'Sim_float_rho'] )

#loop through each of the distance cutoffs
estimate_df_x = []
for curr_x in DIST_TIME_MAX:
    print("X is: " + str(curr_x))
    curr_x_stat = all_stat_dfs[all_stat_dfs['Dist_X_Time'] <= curr_x]

    #loop through each rho value
    for curr_rho in curr_x_stat['Sim_Rho'].unique():
        print(curr_rho)
        curr_rho_stat = curr_x_stat[curr_x_stat['Sim_Rho'] == curr_rho]
        curr_float_rho = all_stat_dfs[all_stat_dfs['Sim_Rho'] == curr_rho]['Sim_float_rho'].unique()[0]

        #loop through each iter group
        for curr_iteration in range(1, NUM_RMSE_GROUPS):
            #get the data for the current rho and iteration
            curr_stat_df = curr_rho_stat[curr_rho_stat['iter_group'] == curr_iteration]

            #sample segregating sites while the D' num is below the threshold
            sample_df = []
            sample_size = 0
            sampled_loci = set()
            unsampled_loci = set(list(curr_stat_df['Locus_1'].unique()) + list(curr_stat_df['Locus_2'].unique().tolist()))
            while sample_size < RATIOS_PER_GROUP and len(unsampled_loci) > 0:
                #sample another locus
                my_sample = np.random.choice(list(unsampled_loci))
                sampled_loci.add(my_sample)
                unsampled_loci.remove(my_sample)
                sample_df = curr_stat_df[curr_stat_df['Locus_1'].isin(sampled_loci) & (curr_stat_df['Locus_2'].isin(sampled_loci))]
                sample_size = len(sample_df)

            curr_stat_df= sample_df

            #get the estimate and fit for the current dataset and sample size
            x_vals = curr_stat_df['Dist_X_Time'].unique()
            lower_fit, upper_fit, estimate_df = plne.bootstrap_rho(curr_stat_df,
                                                                NUM_BOOTSTRAPS)

            lower_conf = np.quantile(estimate_df['Estimated_Rho'], 0.025)
            mid_est = np.quantile(estimate_df['Estimated_Rho'], 0.5)
            upper_conf = np.quantile(estimate_df['Estimated_Rho'], 0.975)

            estimate_df_x.append([
                mid_est, curr_rho, curr_iteration, 
                curr_x, curr_float_rho])  

estimate_df_x = pd.DataFrame(estimate_df_x, columns=["Est_Rho", 
                        'Sim_Rho', 'iter_name', 'x_threshold', 'Sim_float_rho'] )


###################### Getting things set to plot ##############################
params = {'font.size': 8, 'axes.labelsize': 8,'axes.titlesize':8,  'legend.fontsize': 8, 'xtick.labelsize': 8, 'ytick.labelsize': 8}
plt.rcParams.update(params)
linewidth = 1.5

##################### Plotting the MSE for each sample size ###################
#plot the estimates to show how accurate they are
fig, axes = plt.subplots(1, 3, figsize=(7, 2.5))

#First we need to group by the sample size and rho
#Then we can calculate the mean squared error
grouped_ests = estimate_df_seg.groupby(['sample_size', 'Sim_float_rho'])
group_MSE = []
for name,group in grouped_ests:
    truth = [name[1] for x in range(len(group))]
    mse = np.sqrt(mean_squared_error(group['Est_Rho'], truth))/np.mean(group['Est_Rho'])
    string_rho = group['Sim_Rho'].unique()[0]
    group_MSE.append([name[0], name[1], mse, string_rho])

group_MSE = pd.DataFrame(group_MSE, 
    columns=['sample_size', 'Sim_float_rho', 'NRMSE', 'Sim_Rho'])

sns.lineplot(x = 'Sim_float_rho', y = 'NRMSE', 
    data = group_MSE, hue = 'sample_size', ax = axes[2],
    palette=sns.color_palette("icefire", n_colors=len(SEG_LOCI_NUMS)))
axes[2].set_xscale('log')
axes[2].set_ylabel('Normalized RMSE')
axes[2].set_xlabel(r'Simulated Recombination Rate ($\rho$)')
axes[2].legend(title = '# of Loci', loc = 'upper left')

##################### Plotting the MSE for each sample size ###################
#First we need to group by the sample size and rho
#Then we can calculate the mean squared error
grouped_ests = estimate_df_size.groupby(['sample_size', 'Sim_float_rho'])
group_MSE = []
for name,group in grouped_ests:
    truth = [name[1] for x in range(len(group))]
    mse = np.sqrt(mean_squared_error(group['Est_Rho'], truth))/np.mean(group['Est_Rho'])
    string_rho = group['Sim_Rho'].unique()[0]
    group_MSE.append([name[0], name[1], mse, string_rho])

group_MSE = pd.DataFrame(group_MSE, 
    columns=['sample_size', 'Sim_float_rho', 'NRMSE', 'Sim_Rho'])
# print(group)
sns.lineplot(x = 'Sim_float_rho', y = 'NRMSE', 
    data = group_MSE, hue = 'sample_size', ax = axes[1],
    palette=sns.color_palette("icefire", n_colors=len(D_PRIME_NUMS)))
axes[1].set_xscale('log')
axes[1].set_ylabel('Normalized RMSE')
axes[1].set_xlabel(r'Simulated Recombination Rate ($\rho$)')
axes[1].legend(title = '# of LDMs', loc = 'upper left')


##################### Plotting the MSE for each threshold ###################
#First we need to group by the sample size and rho
#Then we can calculate the mean squared error
grouped_ests = estimate_df_x.groupby(['x_threshold', 'Sim_float_rho'])
group_MSE = []
for name,group in grouped_ests:
    truth = [name[1] for x in range(len(group))]
    mse = np.sqrt(mean_squared_error(group['Est_Rho'], truth))/np.mean(group['Est_Rho'])
    string_rho = group['Sim_Rho'].unique()[0]
    group_MSE.append([name[0], name[1], mse, string_rho])

group_MSE = pd.DataFrame(group_MSE, 
    columns=['sample_size', 'Sim_float_rho', 'NRMSE', 'Sim_Rho'])

sns.lineplot(x = 'Sim_float_rho', y = 'NRMSE', 
    data = group_MSE, hue = 'sample_size', ax = axes[0],
   palette=sns.color_palette("icefire", n_colors=len(DIST_TIME_MAX)))
axes[0].set_xscale('log')
axes[0].set_ylabel('Normalized RMSE')
axes[0].set_xlabel(r'Simulated Recombination Rate ($\rho$)')
axes[0].legend(title = r'd$\Delta$t Threshold', loc = 'upper right')
plt.tight_layout()

plt.savefig(outDir + "check_supp_params_" + str(NUM_RMSE_GROUPS) + outTag, dpi = 300)
plt.close()

