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
D_PRIME_NUMS = [2500, 5000, 10000, 25000, 50000]
DIST_TIME_MAX = [5000, 10000, 25000, 50000, 100000, 200000]
NUM_REPS = 200
NUM_GROUPS = 100

#For running on Cluster
dataDir = "/net/feder/vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_10_03_neutral/"
outDir = "/net/feder/vol1/home/evromero/2021_hiv-rec/results/paper/supp_params/"

# dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_10_03_neutral/'
# outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/paper/supp_params/'

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

#Randomly divide the reps into 10 groups
rep_groups = np.array(range(0, NUM_REPS+1))
np.random.shuffle(rep_groups)
rep_groups = np.array_split(rep_groups, NUM_GROUPS)


#Make a dictionary to label each group
group_dict = {}
for i in range(len(rep_groups)):
    for j in rep_groups[i]:
        group_dict[j] = i

#Add the group labels to the dataframe      
group_labels = [group_dict[x] for x in all_stat_dfs['rep']]
all_stat_dfs['iter_group'] = group_labels


#loop through each of the sample sizes
estimate_df_size = []
for curr_size in D_PRIME_NUMS:
    # print("Size is: " + str(curr_size))
    #loop through each rho value
    for curr_rho in all_stat_dfs['Sim_Rho'].unique():
        curr_rho_stat = all_stat_dfs[all_stat_dfs['Sim_Rho'] == curr_rho]
        
        #Add the distance example cutoff
        curr_rho_stat = curr_rho_stat[curr_rho_stat['Dist_X_Time'] <= DIST_EXAMPLE] 

        #loop through each iter group
        for curr_iteration in range(1, NUM_GROUPS+1):
            #get the data for the current rho and iteration
            curr_stat_df = curr_rho_stat[curr_rho_stat['rep'] == curr_iteration]
            
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
            

            #get the estimate and fit for the current dataset and sample size
            x_vals = curr_stat_df['Dist_X_Time'].unique()
            #Get the current estimate
            lower_fit, upper_fit, estimate_df = plne.bootstrap_rho(curr_stat_df,
                                                                NUM_BOOTSTRAPS)
            #04/07/23: making the fit use bootstrapping instead of this single fit
            # coeffs, fit_dat = optimize.curve_fit(plne.neher_leitner, 
            #     curr_stat_df['Dist_X_Time'], curr_stat_df['d_ratio'],
            #     p0 = [0, 0.26, .0000439], maxfev = 10000)

            lower_conf = np.quantile(estimate_df['Estimated_Rho'], 0.025)
            mid_est = np.quantile(estimate_df['Estimated_Rho'], 0.5)
            upper_conf = np.quantile(estimate_df['Estimated_Rho'], 0.975)

            estimate_df_size.append([
                mid_est, curr_data, curr_rho, curr_iteration, 
                curr_data, curr_size])  
            
            # if mid_est < 0:
            #     print("Chaos")
            #     print(curr_rho)
            #     print(curr_size)
            #     binned_rat, binedges, bin_nums = binned_statistic(
            #     curr_stat_df['Dist_X_Time'].to_numpy(), 
            #     curr_stat_df['d_ratio'].to_numpy(), bins = 100)
    
            #     curr_rho_bins = pd.DataFrame({'Dist_X_Time': binedges[1:], 'D_Ratio': binned_rat})
            #     sns.lineplot(data=curr_rho_bins, x='Dist_X_Time', y='D_Ratio')
            #     plt.savefig(outDir + 'Chaos.png')
            #     quit()

estimate_df_size = pd.DataFrame(estimate_df_size, columns=["Est_Rho", 'Dataset', 
                    'Sim_Rho', 'iter_name' , 'data', 'sample_size'] )



#loop through each of the distance cutoffs
estimate_df_x = []
for curr_x in DIST_TIME_MAX:
    # print("X is: " + str(curr_x))
    curr_x_stat = all_stat_dfs[all_stat_dfs['Dist_X_Time'] <= curr_x]

    #loop through each rho value
    for curr_rho in curr_x_stat['Sim_Rho'].unique():
        curr_rho_stat = curr_x_stat[curr_x_stat['Sim_Rho'] == curr_rho]

        #loop through each iter group
        for curr_iteration in range(1, NUM_GROUPS+1):
            #get the data for the current rho and iteration
            curr_stat_df = curr_rho_stat[curr_rho_stat['rep'] == curr_iteration]

            #get the estimate and fit for the current dataset and sample size
            x_vals = curr_stat_df['Dist_X_Time'].unique()
            lower_fit, upper_fit, estimate_df = plne.bootstrap_rho(curr_stat_df,
                                                                NUM_BOOTSTRAPS)


            lower_conf = np.quantile(estimate_df['Estimated_Rho'], 0.025)
            mid_est = np.quantile(estimate_df['Estimated_Rho'], 0.5)
            upper_conf = np.quantile(estimate_df['Estimated_Rho'], 0.975)

            estimate_df_x.append([
                mid_est, curr_data, curr_rho, curr_iteration, 
                curr_data, curr_size, curr_x])  

estimate_df_x = pd.DataFrame(estimate_df_x, columns=["Est_Rho", 'Dataset', 
                        'Sim_Rho', 'iter_name' , 'data', 'sample_size', 'x_threshold'] )


############################# Plotting Estimate Accuracy ######################
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
for entry in estimate_df_x['Sim_Rho']:
    intRhoList.append(rhoDict[entry])
    newStringRho.append(rho_dict_fix_strings[entry])
estimate_df_x['Sim_int_rho'] = intRhoList
estimate_df_x['Sim_Rho'] = newStringRho

#redo the labeling on the rho values from what was used in the simulation names
intRhoList = []
newStringRho = []
for entry in estimate_df_size['Sim_Rho']:
    intRhoList.append(rhoDict[entry])
    newStringRho.append(rho_dict_fix_strings[entry])
estimate_df_size['Sim_int_rho'] = intRhoList
estimate_df_size['Sim_Rho'] = newStringRho


rcParams['mathtext.fontset'] = 'custom'
rcParams['mathtext.rm'] = 'DejaVu Sans'
rcParams['mathtext.it'] = 'DejaVu Sans:italic'

estimate_df = estimate_df_x[estimate_df_x['x_threshold'] == DIST_EXAMPLE]

#plot the estimates to show how accurate they are
sns.set(rc={'figure.figsize':(25,25)}, font_scale = 3, font = '')
sns.set_palette("tab10")
sns.set_style("white")
linewidth = 2

fig, axes = plt.subplots(2)

##################### Plotting the MSE for each sample size ###################
#First we need to group by the sample size and rho
#Then we can calculate the mean squared error
grouped_ests = estimate_df_size.groupby(['sample_size', 'Sim_int_rho'])
group_MSE = []
for name, group in grouped_ests:
    truth = [name[1] for x in range(len(group))]
    mse = np.sqrt(abs(mean_squared_error(group['Est_Rho'], truth)))/abs(np.mean(group['Sim_int_rho']))
    if mse < 0:
        print("The RMSE is")
        print(mse)
        print("The mse is")
        print(mean_squared_error(group['Est_Rho'], truth))
        print("The mean is")
        print(np.mean(group['Est_Rho']) )
    string_rho = group['Sim_Rho'].unique()[0]
    group_MSE.append([name[0], name[1], mse, string_rho])

group_MSE = pd.DataFrame(group_MSE, 
    columns=['sample_size', 'Sim_int_rho', 'NRMSE', 'Sim_Rho'])

sns.lineplot(x = 'Sim_int_rho', y = 'NRMSE', linewidth = linewidth,
    data = group_MSE, hue = 'sample_size', ax = axes[0],
   palette=sns.color_palette("icefire", n_colors=len(D_PRIME_NUMS)))
axes[0].set_xscale('log')
axes[0].set_ylabel('Normalized RMSE')
axes[0].set_xlabel(r'Simulation Value of $\rho$')
axes[0].legend(title = '# of D\' Ratios')


##################### Plotting the MSE for each threshold ###################
#First we need to group by the sample size and rho
#Then we can calculate the mean squared error
grouped_ests = estimate_df_x.groupby(['x_threshold', 'Sim_int_rho'])
group_MSE = []
for name,group in grouped_ests:
    truth = [name[1] for x in range(len(group))]
    mse = np.sqrt(abs(mean_squared_error(group['Est_Rho'], truth)))/abs(np.mean(group['Sim_int_rho']))
    string_rho = group['Sim_Rho'].unique()[0]
    group_MSE.append([name[0], name[1], mse, string_rho])

group_MSE = pd.DataFrame(group_MSE, 
    columns=['sample_size', 'Sim_int_rho', 'NRMSE', 'Sim_Rho'])

sns.lineplot(x = 'Sim_int_rho', y = 'NRMSE', 
    data = group_MSE, hue = 'sample_size', ax = axes[1], 
    linewidth = linewidth, 
    palette=sns.color_palette("icefire", n_colors=len(DIST_TIME_MAX)))
axes[1].set_xscale('log')
axes[1].set_ylabel('Normalized RMSE')
axes[1].set_xlabel(r'Simulation Value of $\rho$')
plt.legend(title = r'd$\Delta$t Threshold')
plt.tight_layout()

plt.savefig(outDir + "supp_params_" + str(NUM_GROUPS) +"_" + str(NUM_BOOTSTRAPS) + ".jpg", dpi = 300)
plt.close()
    