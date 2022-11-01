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
from sklearn.metrics import mean_squared_error
from matplotlib import rcParams
from scipy.stats import linregress

THRESHOLD = 0.2
DIST_EXAMPLE = 200000
DIST_TIME_MAX = [10000, 25000, 50000, 100000, 150000, 200000]
DIST_TIME_MAX = [50000, 200000]
D_PRIME_NUMS = [2500, 5000, 10000, 15000, 20000, 25000]
D_PRIME_NUMS = [25000]
NUM_REPS = 10
NUM_GROUPS = 9
REG_WIN_SIZE = 5000

dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_09_16_neutral_long/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/paper/fig2/'

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
all_stat_dfs['iter_group'] = group_labels


#loop through each of the sample sizes
estimate_df_size = []
for curr_size in D_PRIME_NUMS:
    print("Size is: " + str(curr_size))
    #loop through each rho value
    for curr_rho in all_stat_dfs['Sim_Rho'].unique():
        curr_rho_stat = all_stat_dfs[all_stat_dfs['Sim_Rho'] == curr_rho]



        #loop through each iter group
        for curr_iteration in range(1, NUM_GROUPS+1):
            #get the data for the current rho and iteration
            curr_stat_df = curr_rho_stat[curr_rho_stat['rep'] == curr_iteration]
            
            #sample segregating sites while the D' num is below the threshold
            sample_df = []
            sample_size = 0
            sampled_loci = set()
            unsampled_loci = set(list(curr_stat_df['Locus_1'].unique()) + list(curr_stat_df['Locus_2'].unique().tolist()))
            while sample_size < curr_size:
                #sample another locus
                sampled_loci.add(np.random.choice(list(unsampled_loci)))
                sample_df = curr_stat_df[curr_stat_df['Locus_1'].isin(sampled_loci) & (curr_stat_df['Locus_2'].isin(sampled_loci))]
                sample_size = len(sample_df)
            
            curr_stat_df= sample_df
            x_vals = curr_stat_df['Dist_X_Time'].unique()
            
            num_tests = len(range(int(min(x_vals)), int(max(x_vals)), REG_WIN_SIZE))
            fit_stop = DIST_EXAMPLE

            for i in range(int(min(x_vals)), int(max(x_vals)), REG_WIN_SIZE):
                linreg_df = curr_stat_df[curr_stat_df['Dist_X_Time'].between(i, i+REG_WIN_SIZE)]
                if len(linreg_df) < 10:
                    continue
                mylinreg_results = linregress(linreg_df['Dist_X_Time'], linreg_df['d_ratio'])

                #test to see if the slope is not significantly different from 0
                if (1-mylinreg_results.pvalue) < 0.05:
                    print("stop fitting at " + str(i))
                    print(curr_iteration)
                    print(curr_rho)
                    fit_stop = i
                    break

            #get the estimate and fit for the current dataset and sample size
            curr_stat_df = curr_stat_df[curr_stat_df['Dist_X_Time'] <= fit_stop]
            
         
            coeffs, fit_dat = optimize.curve_fit(plne.neher_leitner, 
                curr_stat_df['Dist_X_Time'], curr_stat_df['d_ratio'],
                p0 = [0, 0.26, .0000439], maxfev = 10000)
            fit_vals = [plne.neher_leitner(x, coeffs[0], coeffs[1], coeffs[2])
                        for x in x_vals]

            #Bin the d' ratios so they are easier to view on the plots
            binned_rat, binedges, bin_nums = binned_statistic(
                curr_stat_df['Dist_X_Time'].to_numpy(), 
                curr_stat_df['d_ratio'].to_numpy(), bins = 100)
            
            my_data = pd.DataFrame(list(zip(binned_rat, binedges[1:])), columns = ['d_ratio', 'Dist_X_Time'])
            my_data.dropna(inplace = True)

            estimate_df_size.append([coeffs[0], coeffs[1], coeffs[2], 
                coeffs[1] * coeffs[2], curr_data, curr_rho, curr_iteration, 
                curr_data, curr_size])  

estimate_df_size = pd.DataFrame(estimate_df_size, columns=["C0", "C1", "C2",
                    "Est_Rho", 'Dataset', 'Sim_Rho', 'iter_name' , 'data', 'sample_size'] )

#loop through each of the distance cutoffs
estimate_df_x = []
for curr_x in DIST_TIME_MAX:
    print("X is: " + str(curr_x))
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

            num_tests = len(range(int(min(x_vals)), int(max(x_vals)), REG_WIN_SIZE))
            fit_stop = curr_x

            for i in range(int(min(x_vals)), int(max(x_vals)), REG_WIN_SIZE):
                linreg_df = curr_stat_df[curr_stat_df['Dist_X_Time'].between(i, i+REG_WIN_SIZE)]
                if len(linreg_df) < 10:
                    continue
                mylinreg_results = linregress(linreg_df['Dist_X_Time'], linreg_df['d_ratio'])

                #test to see if the slope is not significantly different from 0
                if (1-mylinreg_results.pvalue) < 0.05:
                    print("stop fitting at " + str(i))
                    print(curr_iteration)
                    print(curr_rho)
                    fit_stop = i
                    break

            #get the estimate and fit for the current dataset and sample size
            curr_stat_df = curr_stat_df[curr_stat_df['Dist_X_Time'] <= fit_stop]
            print(max(curr_stat_df['Dist_X_Time']))

            coeffs, fit_dat = optimize.curve_fit(plne.neher_leitner, 
                curr_stat_df['Dist_X_Time'], curr_stat_df['d_ratio'],
                p0 = [0, 0.26, .0000439], maxfev = 10000)
            fit_vals = [plne.neher_leitner(x, coeffs[0], coeffs[1], coeffs[2])
                        for x in x_vals]

            #Bin the d' ratios so they are easier to view on the plots
            binned_rat, binedges, bin_nums = binned_statistic(
                curr_stat_df['Dist_X_Time'].to_numpy(), 
                curr_stat_df['d_ratio'].to_numpy(), bins = 100)
            
            my_data = pd.DataFrame(list(zip(binned_rat, binedges[1:])), columns = ['d_ratio', 'Dist_X_Time'])
            my_data.dropna(inplace = True)

            sns.lineplot(x = 'Dist_X_Time', y = 'd_ratio', data = my_data, color = 'black')
            sns.lineplot(x = x_vals, y = fit_vals, color = 'red')
            plt.savefig(outDir + str(curr_rho) + '/fit_results_' + str(curr_iteration) + "_" + str(curr_x) + '.png')
            plt.close()

            estimate_df_x.append([coeffs[0], coeffs[1], coeffs[2], 
                coeffs[1] * coeffs[2], curr_data, curr_rho, curr_iteration, 
                curr_data, len(curr_stat_df), curr_x])  

estimate_df_x = pd.DataFrame(estimate_df_x, columns=["C0", "C1", "C2",
                    "Est_Rho", 'Dataset', 'Sim_Rho', 'iter_name' , 'data', 'sample_size', 'x_threshold'] )



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
sns.set(rc={'figure.figsize':(30,10)}, font_scale = 2, font = '')
fig, axes = plt.subplots(1, 3)
sns.stripplot(x = 'Sim_Rho', y = 'Est_Rho', data = estimate_df, 
    jitter = True, color = 'k', s = 8, ax = axes[0],
    order = [r"$2\times10^{-6}$", r"$10^{-5}$", r"$2\times10^{-5}$", r"$10^{-4}$", r"$2\times10^{-4}$", r"$10^{-3}$"])

# distance across the "X" or "Y" stipplot column to span, in this case 40%
label_width = 0.4

            
for tick, text in zip(axes[0].get_xticks(), axes[0].get_xticklabels()):
    sample_name = text.get_text()  # "X" or "Y"

    #get the float value of rho corresponding with the tick
    rho_val = estimate_df[estimate_df['Sim_Rho'] == sample_name]
    rho_val = rho_val['Sim_int_rho'].unique()[0]

    # plot horizontal lines across the column, centered on the tick
    axes[0].plot([tick-label_width/2, tick+label_width/2], [rho_val, rho_val],
            lw=2, color='k')


axes[0].set_xlabel(r'Simulation Value of $\rho$')
axes[0].set_ylabel(r'Estimated Value of $\rho$')
axes[0].set_ylim(0.000001, 0.01)
axes[0].set_yscale('log')
# axes[0].set_ylim(0.000001, 0.002)


##################### Plotting the MSE for each sample size ###################
#First we need to group by the sample size and rho
#Then we can calculate the mean squared error
grouped_ests = estimate_df_size.groupby(['sample_size', 'Sim_int_rho'])
group_MSE = []
for name,group in grouped_ests:
    truth = [name[1] for x in range(len(group))]
    mse = np.sqrt(mean_squared_error(group['Est_Rho'], truth))/np.mean(group['Est_Rho'])
    string_rho = group['Sim_Rho'].unique()[0]
    group_MSE.append([name[0], name[1], mse, string_rho])

group_MSE = pd.DataFrame(group_MSE, 
    columns=['sample_size', 'Sim_int_rho', 'NRMSE', 'Sim_Rho'])
# print(group)
sns.lineplot(x = 'Sim_int_rho', y = 'NRMSE', 
    data = group_MSE, hue = 'sample_size', ax = axes[1],
   palette=sns.color_palette("icefire", n_colors=len(D_PRIME_NUMS)))
axes[1].set_xscale('log')
axes[1].set_ylabel('Normalized RMSE')
axes[1].set_xlabel(r'Simulation Value of $\rho$')
axes[1].legend(title = '# of D\' Ratios')


##################### Plotting the MSE for each threshold ###################
#First we need to group by the sample size and rho
#Then we can calculate the mean squared error
grouped_ests = estimate_df_x.groupby(['x_threshold', 'Sim_int_rho'])
group_MSE = []
for name,group in grouped_ests:
    truth = [name[1] for x in range(len(group))]
    mse = np.sqrt(mean_squared_error(group['Est_Rho'], truth))/np.mean(group['Est_Rho'])
    string_rho = group['Sim_Rho'].unique()[0]
    group_MSE.append([name[0], name[1], mse, string_rho])

group_MSE = pd.DataFrame(group_MSE, 
    columns=['sample_size', 'Sim_int_rho', 'NRMSE', 'Sim_Rho'])

sns.lineplot(x = 'Sim_int_rho', y = 'NRMSE', 
    data = group_MSE, hue = 'sample_size', ax = axes[2],
   palette=sns.color_palette("icefire", n_colors=len(DIST_TIME_MAX)))
axes[2].set_xscale('log')
axes[2].set_ylabel('Normalized RMSE')
axes[2].set_xlabel(r'Simulation Value of $\rho$')
plt.legend(title = r'd$\Delta$t Threshold')
plt.tight_layout()
    

plt.savefig(outDir + "figure_2_longread.jpg")
plt.close()
    