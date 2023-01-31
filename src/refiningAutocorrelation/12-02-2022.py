import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic
from matplotlib import rcParams
import plot_neher as plne
from scipy import optimize
from scipy.stats import linregress

THRESHOLD = 0.2
DIST_TIME_MAX = 50000
NUM_REPS = 60
NUM_GROUPS = 30
NUM_BOOTSTRAPS = 10
NUM_BINS = 300

#Today I am going to calculate R^2 for each of the fits and maybe fit a line when R^2 is bad

dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_10_03_neutral/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_10_03_neutral/12-02-2022/'

# dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_10_03_neutral/'
# outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_10_03_neutral/'

#First, we will get all of the data and divide it into groups
all_stat_dfs = []

#loop through each of the dataframes for the separate simulations

rho_stat_dfs = []
for curr_data in os.listdir(dataDir):
    #only get the data directories, not hidden files
    if curr_data[0] == '.':
        continue

    #get the information for the current run
    run_info = curr_data.split('_')
    sim_rho = run_info[1]
    sim_rho = sim_rho[3:]
    rep = run_info[-1]
    if int(rep[3:]) >= NUM_REPS:
        continue

    #get the dataframe for the current run
    d_ratio_file = dataDir + curr_data + "/linkage/d_ratio"
    stat_df = pd.read_pickle(d_ratio_file)
    stat_df['rep'] = int(rep[3:])
    stat_df['Sim_Rho'] = sim_rho
    rho_stat_dfs.append(stat_df)
rho_stat_dfs = pd.concat(rho_stat_dfs)

#Randomly divide the reps into 10 groups
rep_groups = np.array(range(0, NUM_REPS))
np.random.shuffle(rep_groups)
rep_groups = np.array_split(rep_groups, NUM_GROUPS)


#Make a dictionary to label each group
group_dict = {}
for i in range(len(rep_groups)):
    for j in rep_groups[i]:
        group_dict[j] = i
        
group_labels = [group_dict[x] for x in rho_stat_dfs['rep']]
rho_stat_dfs['iter_group'] = group_labels
all_stat_dfs.append(rho_stat_dfs)

all_stat_dfs = pd.concat(all_stat_dfs, ignore_index=True)
all_stat_dfs = all_stat_dfs[all_stat_dfs['Dist_X_Time'] <= DIST_TIME_MAX]
all_stat_dfs = all_stat_dfs[all_stat_dfs['d_ratio'] >= THRESHOLD]

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
float_to_str = {0.001 : r"$10^{-3}$",
            0.0001 : r"$10^{-4}$",
            0.0002 : r"$2\times10^{-4}$",
            0.00001 : r"$10^{-5}$",
            0.00002 : r"$2\times10^{-5}$",
            0.000002 : r"$2\times10^{-6}$"}

#redo the labeling on the rho values from what was used in the simulation names
intRhoList = []
newStringRho = []
for entry in all_stat_dfs['Sim_Rho']:
    intRhoList.append(rhoDict[entry])
    newStringRho.append(rho_dict_fix_strings[entry])
all_stat_dfs['Sim_float_rho'] = intRhoList
all_stat_dfs['Sim_Rho'] = newStringRho
all_stat_dfs['dist'] = all_stat_dfs['Locus_2'] - all_stat_dfs['Locus_1']
print(all_stat_dfs)

####################### Plotting the Ddelta t distribution ####################
sns.histplot(data = all_stat_dfs, stat = 'density', x = 'Dist_X_Time', color = 'blue', bins = NUM_BINS)

plt.savefig(outDir + "Ddelta_t_dist.png", dpi = 300)
plt.close()

sns.scatterplot(data = all_stat_dfs, x = 'dist', y = 'Time_Diff')
plt.ylim(0,700)

plt.savefig(outDir + "dist_and_time.png", dpi = 300)
plt.close()

####### Now we can check the estimates with the downsampled distribution ######

########################## Estimating recombination rates #####################
#loop through each of the distance cutoffs
all_stat_dfs = all_stat_dfs[all_stat_dfs['Dist_X_Time'] <= DIST_TIME_MAX]
all_estimate_df = []
all_best_fits = []

#loop through each rho value
for curr_rho in all_stat_dfs['Sim_Rho'].unique():
    curr_rho_stat = all_stat_dfs[all_stat_dfs['Sim_Rho'] == curr_rho]
    curr_float_rho = curr_rho_stat['Sim_float_rho'].unique()[0]

    #loop through each iter group
    for curr_iteration in range(0, NUM_GROUPS):
        #get the data for the current rho and iteration
        curr_stat_df = curr_rho_stat[curr_rho_stat['iter_group'] == curr_iteration]

        #Now downsample the group to the coverage of the lowest bin
        #Bin the d' ratios so they are easier to view on the plots
        bin_count, binedges, bin_nums = binned_statistic(
        curr_stat_df['Dist_X_Time'].to_numpy(), 
        curr_stat_df['d_ratio'].to_numpy(), bins = 100, statistic= 'count')


        #get the estimate and fit for the current dataset and sample size
        x_vals = curr_stat_df['Dist_X_Time'].unique()
        coeffs, fit_dat = optimize.curve_fit(plne.neher_leitner, 
            curr_stat_df['Dist_X_Time'], curr_stat_df['d_ratio'],
            p0 = [0, 0.26, .0000439], maxfev = 100000)
        fit_vals = [plne.neher_leitner(x, coeffs[0], coeffs[1], coeffs[2])
                    for x in x_vals]
        
        line_fit = linregress(curr_stat_df['Dist_X_Time'], curr_stat_df['d_ratio'])
        short_dist = curr_stat_df[curr_stat_df['Dist_X_Time'] <= 2000]
        short_dist_fit = linregress(short_dist['Dist_X_Time'], short_dist['d_ratio'])
        # print("Short dist fit: " + str(short_dist_fit.rvalue**2))
        # print("Full dist fit: " + str(line_fit.rvalue**2))



        #calculate the r2 value for the fit
        residuals, res_bins, res_bin_nums = binned_statistic(curr_stat_df['Dist_X_Time'].to_numpy(),
            curr_stat_df['d_ratio'] - [plne.neher_leitner(x, coeffs[0], coeffs[1], coeffs[2])
                    for x in curr_stat_df['Dist_X_Time']], bins = 100, statistic = 'mean')
        tot_ss, tot_bins, tot_bin_nums = binned_statistic(curr_stat_df['Dist_X_Time'].to_numpy(), curr_stat_df['d_ratio'].to_numpy(), bins = 100, statistic= 'mean')
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((tot_ss - np.mean(curr_stat_df['d_ratio']))**2)
        r2 = 1 - (ss_res / ss_tot)   

        #choose the fit with the best r2 value
        which_fit = np.argmax([short_dist_fit.rvalue**2, line_fit.rvalue**2, r2])

        best_fit_df = []
        if which_fit == 0:
            print("Short dist")
            best_fit_df.append([short_dist_fit.slope, short_dist_fit.rvalue**2])
        elif which_fit == 1:
            print("line")
            best_fit_df.append([line_fit.slope, line_fit.rvalue**2])
        else:
            best_fit_df.append([coeffs[1]*coeffs[2], r2])

        best_fit_df = pd.DataFrame(best_fit_df, columns = ['Estimated_Rho', 'Fit_R2'])
        best_fit_df['Sim_Rho'] = curr_rho
        all_best_fits.append(best_fit_df)
        print(which_fit)
        print('************************************************')

        #Bin the d' ratios so they are easier to view on the plots
        binned_rat, binedges, bin_nums = binned_statistic(
            curr_stat_df['Dist_X_Time'].to_numpy(), 
            curr_stat_df['d_ratio'].to_numpy(), bins = 100)

        sns.lineplot(y = binned_rat, x = binedges[1:], color = 'black')
        sns.lineplot(x = x_vals, y = fit_vals, color = 'red')
        plt.savefig(outDir + 'fits/fit_results_r2_' + str(curr_float_rho) + "_" + str(curr_iteration) + '.png')
        plt.close()

        estimate_df = pd.DataFrame([[coeffs[0], coeffs[1], coeffs[2], coeffs[1] * coeffs[2], r2]], columns = ['c0', 'c1', 'c2', 'Estimated_Rho', 'Fit_R2'])

        estimate_df['Group'] = curr_iteration
        estimate_df['Sim_Rho'] = curr_rho
        all_estimate_df.append(estimate_df)


all_estimate_df = pd.concat(all_estimate_df, ignore_index=True)
all_best_fits = pd.concat(all_best_fits, ignore_index=True)


#To get the comparison, we need to convert the rho values to floats
rhoDict = {r"$10^{-3}$" : 0.001,
            r"$10^{-4}$" : 0.0001,
            r"$2\times10^{-4}$" : 0.0002,
            r"$10^{-5}$" : 0.00001,
            r"$2\times10^{-5}$" : 0.00002,
            r"$2\times10^{-6}$" : 0.000002}


#redo the labeling on the rho values from what was used in the simulation names
intRhoList = []
newStringRho = []
for entry in all_estimate_df['Sim_Rho']:
    intRhoList.append(rhoDict[entry])
all_estimate_df['Sim_float_rho'] = intRhoList

#redo the labeling on the rho values from what was used in the simulation names
intRhoList = []
newStringRho = []
for entry in all_best_fits['Sim_Rho']:
    intRhoList.append(rhoDict[entry])
all_best_fits['Sim_float_rho'] = intRhoList



########################## Plotting the accuracy of the estimates #############
#plot the estimates to show how accurate they are
sns.set(rc={'figure.figsize':(10,10)}, font_scale = 2, font = '')
sns.set_palette("tab10")
sns.set_style("white")

fig = sns.stripplot(x = 'Sim_Rho', y = 'Estimated_Rho', data = all_estimate_df, 
    jitter = True, hue = 'Fit_R2', s = 8, alpha = 1,
    order = [r"$2\times10^{-6}$", r"$10^{-5}$", r"$2\times10^{-5}$", r"$10^{-4}$", r"$2\times10^{-4}$" , r"$10^{-3}$"])

# distance across the "X" or "Y" stipplot column to span, in this case 40%
label_width = 0.4

#need to save the fig to force the axis to be drawn
plt.savefig(outDir + "accuracy_goodness_of_fit.png", dpi = 300)
            
for tick, text in zip(fig.get_xticks(), fig.get_xticklabels()):
    sample_name = text.get_text()  # "X" or "Y"

    #get the float value of rho corresponding with the tick
    rho_val = all_estimate_df[all_estimate_df['Sim_Rho'] == sample_name]
    rho_val = rho_val['Sim_float_rho'].unique()[0]

    # plot horizontal lines across the column, centered on the tick
    fig.plot([tick-label_width/2, tick+label_width/2], [rho_val, rho_val],
            lw=2, color='k')

plt.ylabel(r'Estimated Value of $\rho$')
plt.xlabel(r'Simulation Value of $\rho$')
plt.yscale('log')
plt.savefig(outDir + "accuracy_goodness_of_fit.png", dpi = 300)
plt.close()

#plot the estimates to show how accurate they are
sns.set(rc={'figure.figsize':(10,10)}, font_scale = 2, font = '')
sns.set_palette("tab10")
sns.set_style("white")

fig = sns.stripplot(x = 'Sim_Rho', y = 'Estimated_Rho', data = all_best_fits, 
    jitter = True, hue = 'Fit_R2', s = 8, alpha = 1,
    order = [r"$2\times10^{-6}$", r"$10^{-5}$", r"$2\times10^{-5}$", r"$10^{-4}$", r"$2\times10^{-4}$" , r"$10^{-3}$"])

# distance across the "X" or "Y" stipplot column to span, in this case 40%
label_width = 0.4

#need to save the fig to force the axis to be drawn
plt.savefig(outDir + "accuracy_goodness_of_choice.png", dpi = 300)
            
for tick, text in zip(fig.get_xticks(), fig.get_xticklabels()):
    sample_name = text.get_text()  # "X" or "Y"

    #get the float value of rho corresponding with the tick
    rho_val = all_estimate_df[all_estimate_df['Sim_Rho'] == sample_name]
    rho_val = rho_val['Sim_float_rho'].unique()[0]

    # plot horizontal lines across the column, centered on the tick
    fig.plot([tick-label_width/2, tick+label_width/2], [rho_val, rho_val],
            lw=2, color='k')

plt.ylabel(r'Estimated Value of $\rho$')
plt.xlabel(r'Simulation Value of $\rho$')
plt.yscale('log')
plt.savefig(outDir + "accuracy_goodness_of_fit_choice.png", dpi = 300)
plt.close()