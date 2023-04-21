import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import zaniniUtil as zu
from scipy.stats import binned_statistic
import plot_neher as plne
import autocorrelation as autocorr
from scipy import optimize
from sklearn.metrics import mean_squared_error
from matplotlib import rcParams

THRESHOLD = 0.2
DIST_EXAMPLE = 50000
NUM_BOOTSTRAPS = 10#00
D_PRIME_NUMS = [25000]
NUM_REPS = 50#0
NUM_GROUPS = 10#0

#For running on Cluster
# dataDir = "/net/feder/vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_10_03_neutral/"
# outDir = "/net/feder/vol1/home/evromero/2021_hiv-rec/results/paper/supp_params/"

dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_10_03_neutral/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_10_03_neutral/'

#First, we will get all of the data and divide it into groups
all_stat_dfs = []

#loop through each of the dataframes for the separate simulations
for curr_data in os.listdir(dataDir):
    #only get the data directories, not hidden files
    if curr_data[0] == '.':
        continue

    #get the information for the current run
    run_info = curr_data.split('_')
    sim_rho = run_info[1]
    sim_rho = sim_rho[3:]
    rep = run_info[-1]
    if int(rep[3:]) > NUM_REPS:
        continue

    #get the dataframe for the current run
    d_ratio_file = dataDir + curr_data + "/linkage/d_ratio"
    stat_df = pd.read_pickle(d_ratio_file)
    stat_df['rep'] = int(rep[3:])
    stat_df['Sim_Rho'] = sim_rho
    all_stat_dfs.append(stat_df)
all_stat_dfs = pd.concat(all_stat_dfs)

#Filter the groups so the first d' value is greater than 0.2
all_stat_dfs = all_stat_dfs[all_stat_dfs['d_i'] > 0.2]

##Make a histogram of the number of D' ratios in each simulation
ratio_nums = []
for name, group in all_stat_dfs.groupby(['Sim_Rho', 'rep']):
    ratio_nums.append([name[0], name[1], len(group)])

ratio_nums = pd.DataFrame(ratio_nums, columns=['Sim_Rho', 'rep', 'num_ratios'])
sns.histplot(data=ratio_nums, x='num_ratios')
plt.xlabel('Number of D\' Ratios')
plt.ylabel('Number of Simulations')
plt.savefig(outDir + 'num_ratios.png')


#Randomly divide the reps into 10 groups
rep_groups = np.array(range(0, NUM_REPS+1))
np.random.shuffle(rep_groups)
rep_groups = np.array_split(rep_groups, NUM_GROUPS)
print(rep_groups)

#Make a dictionary to label each group
group_dict = {}
for i in range(len(rep_groups)):
    for j in rep_groups[i]:
        group_dict[j] = i

#Add the group labels to the dataframe      
group_labels = [group_dict[x] for x in all_stat_dfs['rep']]
all_stat_dfs['iter_group'] = group_labels
print(group_dict)


#loop through each of the sample sizes
estimate_df_size = []
for curr_size in D_PRIME_NUMS:
    print("Size is: " + str(curr_size))
    #loop through each rho value
    for curr_rho in all_stat_dfs['Sim_Rho'].unique():
        curr_rho_stat = all_stat_dfs[all_stat_dfs['Sim_Rho'] == curr_rho]
        
        #Add the distance example cutoff
        curr_rho_stat = curr_rho_stat[curr_rho_stat['Dist_X_Time'] <= DIST_EXAMPLE] 

        #loop through each iter group
        for curr_iteration in range(1, NUM_GROUPS):
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
            if len(unsampled_loci) == 0:
                print("No more loci to sample")
                print(curr_iteration)
                quit()

            curr_stat_df= sample_df
            print(len(curr_stat_df))

            #get the estimate and fit for the current dataset and sample size
            x_vals = curr_stat_df['Dist_X_Time'].unique()
            #Get the current estimate
            lower_fit, upper_fit, estimate_df = plne.bootstrap_rho(curr_stat_df,
                                                                NUM_BOOTSTRAPS)

            lower_conf = np.quantile(estimate_df['Estimated_Rho'], 0.025)
            mid_est = np.quantile(estimate_df['Estimated_Rho'], 0.5)
            upper_conf = np.quantile(estimate_df['Estimated_Rho'], 0.975)

            estimate_df_size.append([
                mid_est, curr_data, curr_rho, curr_iteration, 
                curr_data, curr_size])  
            
estimate_df_size = pd.DataFrame(estimate_df_size, columns=["Est_Rho", 'Dataset', 
                    'Sim_Rho', 'iter_name' , 'data', 'sample_size'] )

########### Now just use one group
print("Now just use one group")
estimate_df_one = []
for curr_rho in all_stat_dfs['Sim_Rho'].unique():
        curr_rho_stat = all_stat_dfs[all_stat_dfs['Sim_Rho'] == curr_rho]
        
        #Add the distance example cutoff
        curr_rho_stat = curr_rho_stat[curr_rho_stat['Dist_X_Time'] <= DIST_EXAMPLE] 

        #loop through each iter group
        for curr_rep in range(0, 10):
            #get the data for the current rho and iteration
            curr_stat_df = curr_rho_stat[curr_rho_stat['rep'] == curr_rep]
        
            print(len(curr_stat_df))
            #get the estimate and fit for the current dataset and sample size
            x_vals = curr_stat_df['Dist_X_Time'].unique()
            #Get the current estimate
            lower_fit, upper_fit, estimate_df = plne.bootstrap_rho(curr_stat_df,
                                                                NUM_BOOTSTRAPS)

            lower_conf = np.quantile(estimate_df['Estimated_Rho'], 0.025)
            mid_est = np.quantile(estimate_df['Estimated_Rho'], 0.5)
            upper_conf = np.quantile(estimate_df['Estimated_Rho'], 0.975)

            estimate_df_one.append([
                mid_est, curr_data, curr_rho, curr_rep, 
                curr_data, curr_size])  
            
estimate_df_one = pd.DataFrame(estimate_df_size, columns=["Est_Rho", 'Dataset', 
                    'Sim_Rho', 'iter_name' , 'data', 'sample_size'] )

############################### Make anothere figure with 20k and 50k loci accuracy side by side ###############################
# #Plot our estimates against each other 
#make the rho values ints rather than strings
rhoDict = {"0.001" : 0.001,
            "1e-04" : 0.0001,
            "2e-04" : 0.0002,
            "1e-05" : 0.00001,
            "2e-05" : 0.00002,
            "2e-06" : 0.000002}
rho_dict_fix_strings = { "0.001" : r"$10^{-3}$",
                        "1e-04" : r"$10^{-4}$",
                        "2e-04" : r"$2\times10^{-4}$",
                        "1e-05" : r"$10^{-5}$",
                        "2e-05" : r"$2\times10^{-5}$",
                        "2e-06" : r"$2\times10^{-6}$"}

#redo the labeling on the rho values from what was used in the simulation names
intRhoList = []
newStringRho = []
for entry in estimate_df_one['Sim_Rho']:
    intRhoList.append(rhoDict[entry])
    newStringRho.append(rho_dict_fix_strings[entry])
estimate_df_one['Sim_int_rho'] = intRhoList
estimate_df_one['Sim_Rho'] = newStringRho

#redo the labeling on the rho values from what was used in the simulation names
intRhoList = []
newStringRho = []
for entry in estimate_df_size['Sim_Rho']:
    intRhoList.append(rhoDict[entry])
    newStringRho.append(rho_dict_fix_strings[entry])
estimate_df_size['Sim_int_rho'] = intRhoList
estimate_df_size['Sim_Rho'] = newStringRho

#plot the estimates to show how accurate they are
sns.set(rc={'figure.figsize':(20,10)}, font_scale = 2, font = '')
fig, axes = plt.subplots(1, 2)
sns.stripplot(x = 'Sim_Rho', y = 'Est_Rho', data = estimate_df_size, 
    jitter = True, color = 'k', s = 8, ax = axes[0],
    order = [r"$2\times10^{-6}$", r"$10^{-5}$", r"$2\times10^{-5}$", r"$10^{-4}$", r"$2\times10^{-4}$", r"$10^{-3}$"])

# distance across the "X" or "Y" stipplot column to span, in this case 40%
label_width = 0.4

for tick, text in zip(axes[0].get_xticks(), axes[0].get_xticklabels()):
    sample_name = text.get_text()  # "X" or "Y"

    #get the float value of rho corresponding with the tick
    rho_val = estimate_df_size[estimate_df_size['Sim_Rho'] == sample_name]
    rho_val = rho_val['Sim_int_rho'].unique()[0]

    # plot horizontal lines across the column, centered on the tick
    axes[0].plot([tick-label_width/2, tick+label_width/2], [rho_val, rho_val],
            lw=2, color='k')


axes[0].set_xlabel(r'Simulation Value of $\rho$')
axes[0].set_ylabel(r'Estimated Value of $\rho$')
axes[0].set_ylim(0.000001, 0.01)
axes[0].set_yscale('log')
axes[0].set_title('Two Sampled Reps')

sns.stripplot(x = 'Sim_Rho', y = 'Est_Rho', data = estimate_df_one, 
    jitter = True, color = 'k', s = 8, ax = axes[1],
    order = [r"$2\times10^{-6}$", r"$10^{-5}$", r"$2\times10^{-5}$", r"$10^{-4}$", r"$2\times10^{-4}$", r"$10^{-3}$"])

for tick, text in zip(axes[1].get_xticks(), axes[1].get_xticklabels()):
    sample_name = text.get_text()  # "X" or "Y"

    #get the float value of rho corresponding with the tick
    rho_val = estimate_df_one[estimate_df_one['Sim_Rho'] == sample_name]
    rho_val = rho_val['Sim_int_rho'].unique()[0]

    # plot horizontal lines across the column, centered on the tick
    axes[1].plot([tick-label_width/2, tick+label_width/2], [rho_val, rho_val],
            lw=2, color='k')
    
axes[1].set_xlabel(r'Simulation Value of $\rho$')
axes[1].set_ylabel(r'Estimated Value of $\rho$')
axes[1].set_ylim(0.000001, 0.01)
axes[1].set_yscale('log')
axes[1].set_title('One Rep')

plt.tight_layout()

plt.savefig(outDir + "accuracy_comparison_one_two_" + str(NUM_GROUPS) + ".jpg")
plt.close()