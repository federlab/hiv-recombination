import sys
from turtle import filling
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
# from sklearn.metrics import mean_squared_error
from matplotlib import rcParams

THRESHOLD = 0.2
D_PRIME_NUMS = [5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000]
NUM_REPS = 160
NUM_GROUPS = 10

dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_09_07_neutral/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_09_07_neutral/'

 
all_stat_dfs = []
estimate_df = []

#loop through each of the dataframes for the separeate simulations
for curr_data in os.listdir(dataDir):
    #only get the data directories, not hidden files
    if curr_data[0] == '.':
        continue

    #get the information for the current run
    run_info = curr_data.split('_')
    sim_rho = run_info[1]
    sim_rho = sim_rho[3:]
    rep = run_info[-1]

    #get the dataframe for the current run
    d_ratio_file = dataDir + curr_data + "/linkage/d_ratio"
    stat_df = pd.read_pickle(d_ratio_file)
    stat_df['rep'] = int(rep[3:])
    stat_df['Sim_Rho'] = sim_rho
    all_stat_dfs.append(stat_df)
all_stat_dfs = pd.concat(all_stat_dfs)


#Randomly divide the reps into 10 groups
rep_groups = np.array(range(0, NUM_REPS))
np.random.shuffle(rep_groups)
rep_groups = np.array_split(rep_groups, NUM_GROUPS)

#Make a dictionary to label each group
group_dict = {}
for i in range(len(rep_groups)):
    for j in rep_groups[i]:
        group_dict[j] = i
        
group_labels = [group_dict[x] for x in all_stat_dfs['rep']]
print(group_labels)
all_stat_dfs['iter_group'] = group_labels

#loop through each rho value
print(all_stat_dfs['Sim_Rho'].unique())
for curr_rho in all_stat_dfs['Sim_Rho'].unique():
    curr_rho_stat = all_stat_dfs[all_stat_dfs['Sim_Rho'] == curr_rho]
    print(curr_rho)

    #loop through each iter group
    for curr_iteration in range(NUM_GROUPS):
        print(curr_iteration)
        #get the data for the current rho and iteration
        curr_stat_df = curr_rho_stat[curr_rho_stat['iter_group'] == curr_iteration]
        # print(curr_stat_df)
        # print(curr_iteration)

        #get the estimate and fit for the current dataset and sample size
        x_vals = curr_stat_df['Dist_X_Time'].unique()
        coeffs, fit_dat = optimize.curve_fit(plne.neher_leitner, 
            curr_stat_df['Dist_X_Time'], curr_stat_df['d_ratio'],
            p0 = [0, 0.26, .0000439], maxfev = 10000)
        fit_vals = [plne.neher_leitner(x, coeffs[0], coeffs[1], coeffs[2])
                    for x in x_vals]

        #Bin the d' ratios so they are easier to view on the plots
        binned_rat, binedges, bin_nums = binned_statistic(
            curr_stat_df['Dist_X_Time'].to_numpy(), 
            curr_stat_df['d_ratio'].to_numpy(), bins = 100)

        estimate_df.append([coeffs[0], coeffs[1], coeffs[2], 
            coeffs[1] * coeffs[2], curr_data, curr_rho, curr_iteration, 
            curr_data, len(curr_stat_df)])  

estimate_df = pd.DataFrame(estimate_df, columns=["C0", "C1", "C2",
                     "Est_Rho", 'Dataset', 'Sim_Rho', 'iter_name' , 'data', 'sample_size'] )
print(estimate_df)
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
for entry in estimate_df['Sim_Rho']:
    intRhoList.append(rhoDict[entry])
    newStringRho.append(rho_dict_fix_strings[entry])
estimate_df['Sim_int_rho'] = intRhoList
estimate_df['Sim_Rho'] = newStringRho


rcParams['mathtext.fontset'] = 'custom'
rcParams['mathtext.rm'] = 'DejaVu Sans'
rcParams['mathtext.it'] = 'DejaVu Sans:italic'


#plot the estimates to show how accurate they are
sns.set(rc={'figure.figsize':(20,10)}, font_scale = 2, font = '')
fig, axes = plt.subplots(1, 1)
sns.stripplot(x = 'Sim_Rho', y = 'Est_Rho', data = estimate_df, 
    jitter = True, color = 'k', s = 8, ax = axes,
    order = [r"$2\times10^{-6}$", r"$10^{-5}$", r"$2\times10^{-5}$", r"$10^{-4}$", r"$2\times10^{-4}$", r"$10^{-3}$"])

# distance across the "X" or "Y" stipplot column to span, in this case 40%
label_width = 0.4

            
for tick, text in zip(axes.get_xticks(), axes.get_xticklabels()):
    sample_name = text.get_text()  # "X" or "Y"

    #get the float value of rho corresponding with the tick
    rho_val = estimate_df[estimate_df['Sim_Rho'] == sample_name]
    rho_val = rho_val['Sim_int_rho'].unique()[0]

    # plot horizontal lines across the column, centered on the tick
    axes.plot([tick-label_width/2, tick+label_width/2], [rho_val, rho_val],
            lw=2, color='k')


axes.set_xlabel(r'Simulation Value of $\rho$')
axes.set_ylabel(r'Estimated Value of $\rho$')
axes.set_ylim(0.000001, 0.01)
axes.set_yscale('log')

plt.savefig(outDir + "accuracy.jpg")
plt.close()

sns.histplot(data = estimate_df, x = 'sample_size')
plt.savefig(outDir + "sample_size.jpg")

# ##################### Plotting the MSE for each sample size ###################
# #First we need to group by the sample size and rho
# #Then we can calculate the mean squared error
# grouped_ests = estimate_df.groupby(['Sample Size', 'Sim_int_rho'])
# group_MSE = []
# for name,group in grouped_ests:
#     truth = [name[1] for x in range(len(group))]
#     mse = np.sqrt(mean_squared_error(group['Est_Rho'], truth))/np.mean(group['Est_Rho'])
#     string_rho = group['Sim_Rho'].unique()[0]
#     group_MSE.append([name[0], name[1], mse, string_rho])

# group_MSE = pd.DataFrame(group_MSE, 
#     columns=['Sample Size', 'Sim_int_rho', 'NRMSE', 'Sim_Rho'])

# sns.lineplot(x = 'Sim_int_rho', y = 'NRMSE', 
#     data = group_MSE, hue = 'Sample Size', ax = axes[1],
#     palette=sns.color_palette("colorblind", n_colors=4))
# axes[1].set_xscale('log')
# axes[1].set_ylabel('Normalized RMSE')
# axes[1].set_xlabel(r'Simulation Value of $\rho$')
# plt.legend(title = 'Segregating Loci')
# plt.tight_layout()

