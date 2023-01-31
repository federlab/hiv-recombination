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
from scipy.stats import linregress

THRESHOLD = 0.2
DIST_TIME_MAX = 100000
NUM_REPS = 10
NUM_GROUPS = 9
NUM_BOOTSTRAPS = 10
NUM_BINS = 300

#Today I am going to calculate R^2 for each of the fits and maybe fit a line when R^2 is bad

dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_09_16_neutral_long/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_09_16_neutral_long/01-30-2023/'

# dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_10_03_neutral/'
# outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_10_03_neutral/01-30-2023/'

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
all_stat_dfs = all_stat_dfs[all_stat_dfs['d_i'] >= THRESHOLD]

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

########################## Estimating recombination rates #####################
#loop through each of the distance cutoffs
all_stat_dfs = all_stat_dfs[all_stat_dfs['Dist_X_Time'] <= DIST_TIME_MAX]
print(max(all_stat_dfs['Dist_X_Time']))
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

        long_dist_fit = linregress(long_dist['Dist_X_Time'], long_dist['d_ratio'])
        line_fit = linregress(curr_stat_df['Dist_X_Time'], curr_stat_df['d_ratio'])
        mid_dist = curr_stat_df[curr_stat_df['Dist_X_Time'] <= 10000]
        mid_dist_fit = linregress(mid_dist['Dist_X_Time'], mid_dist['d_ratio'])
        short_dist = curr_stat_df[curr_stat_df['Dist_X_Time'] <= 2000]
        short_dist_fit = linregress(short_dist['Dist_X_Time'], short_dist['d_ratio'])
        super_short_dist = curr_stat_df[curr_stat_df['Dist_X_Time'] <= 1000]
        super_short_dist_fit = linregress(super_short_dist['Dist_X_Time'], super_short_dist['d_ratio'])
        shortest_dist = curr_stat_df[curr_stat_df['Dist_X_Time'] <= 500]
        shortest_dist_fit = linregress(shortest_dist['Dist_X_Time'], shortest_dist['d_ratio'])


        #choose the fit with the best r2 value
        which_fit = np.argmax([short_dist_fit.rvalue**2, line_fit.rvalue**2, mid_dist_fit.rvalue**2, super_short_dist_fit.rvalue**2, shortest_dist_fit.rvalue**2, long_dist_fit.rvalue**2])

        best_fit_df = []
        if which_fit == 0:
            print("Short dist")
            fit_vals = short_dist_fit.slope * x_vals + short_dist_fit.intercept
        elif which_fit == 1:
            print("line")
            fit_vals = line_fit.slope * x_vals + line_fit.intercept
        elif which_fit == 2:
            print("mid dist")
            fit_vals = mid_dist_fit.slope * x_vals + mid_dist_fit.intercept
        elif which_fit == 3:
            print("super short dist")
            fit_vals = super_short_dist_fit.slope * x_vals + super_short_dist_fit.intercept
        elif which_fit == 4:
            print("shortest dist")
            fit_vals = shortest_dist_fit.slope * x_vals + shortest_dist_fit.intercept
        elif which_fit == 5:
            print("long dist")
            fit_vals = long_dist_fit.slope * x_vals + long_dist_fit.intercept

        best_fit_df.append([shortest_dist_fit.slope, shortest_dist_fit.rvalue**2, 500])
        best_fit_df.append([super_short_dist_fit.slope, short_dist_fit.rvalue**2, 1000])
        best_fit_df.append([short_dist_fit.slope, short_dist_fit.rvalue**2, 2000])
        best_fit_df.append([line_fit.slope, line_fit.rvalue**2, 50000])
        best_fit_df.append([mid_dist_fit.slope, mid_dist_fit.rvalue**2, 10000])
        best_fit_df.append([long_dist_fit.slope, long_dist_fit.rvalue**2, 100000])


        best_fit_df = pd.DataFrame(best_fit_df, columns = ['Estimated_Rho', 'Fit_R2', 'Fit_Type'])
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
for entry in all_best_fits['Sim_Rho']:
    intRhoList.append(rhoDict[entry])
all_best_fits['Sim_float_rho'] = intRhoList



########################## Plotting the accuracy of the estimates #############
#plot the estimates to show how accurate they are
sns.set(rc={'figure.figsize':(10,10)}, font_scale = 2, font = '')
sns.set_palette("tab10")
sns.set_style("white")

fig = sns.stripplot(x = 'Sim_Rho', y = 'Estimated_Rho', data = all_best_fits, 
    jitter = True, hue = 'Fit_Type', s = 8, alpha = 1,
    order = [r"$2\times10^{-6}$", r"$10^{-5}$", r"$2\times10^{-5}$", r"$10^{-4}$", r"$2\times10^{-4}$" , r"$10^{-3}$"],
    palette=sns.color_palette("coolwarm", n_colors=len(all_best_fits['Fit_Type'].unique())))

# distance across the "X" or "Y" stipplot column to span, in this case 40%
label_width = 0.4

#need to save the fig to force the axis to be drawn
plt.savefig(outDir + "line_fits_both.png", dpi = 300)
            
for tick, text in zip(fig.get_xticks(), fig.get_xticklabels()):
    sample_name = text.get_text()  # "X" or "Y"

    #get the float value of rho corresponding with the tick
    rho_val = all_best_fits[all_best_fits['Sim_Rho'] == sample_name]
    rho_val = rho_val['Sim_float_rho'].unique()[0]

    # plot horizontal lines across the column, centered on the tick
    fig.plot([tick-label_width/2, tick+label_width/2], [rho_val, rho_val],
            lw=2, color='k')

plt.ylabel(r'Estimated Value of $\rho$')
plt.xlabel(r'Simulation Value of $\rho$')
plt.yscale('log')
plt.savefig(outDir + "line_fits_both.png", dpi = 300)
plt.close()