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

DI_THRESH = 0.2
DIST_TIME_MAX = 50000
NUM_REPS = 2
NUM_GROUPS = 1
NUM_BOOTSTRAPS = 10
NUM_BINS = 100

#Today I am trying to resample the ddelta t values to see if that is affecting the fit

dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_10_03_neutral/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_10_03_neutral/01-26-2023/'

# dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_10_03_neutral/'
# outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_10_03_neutral/01-26-2023/'

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
all_stat_dfs = all_stat_dfs[all_stat_dfs['d_i'] > DI_THRESH]
# all_stat_dfs = all_stat_dfs[all_stat_dfs['d_ratio'] >= THRESHOLD]

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
all_fit_df = []

#loop through each rho value
for curr_rho in all_stat_dfs['Sim_Rho'].unique():
    curr_rho_stat = all_stat_dfs[all_stat_dfs['Sim_Rho'] == curr_rho]
    curr_float_rho = curr_rho_stat['Sim_float_rho'].unique()[0]

    #loop through each iter group
    for curr_iteration in range(0, NUM_GROUPS):
        #get the data for the current rho and iteration
        curr_stat_df = curr_rho_stat[curr_rho_stat['iter_group'] == curr_iteration]

        # # #Try binning by unique loci
        # unique_loci = 

        #Now downsample the group to the coverage of the lowest bin
        #Bin the d' ratios so they are easier to view on the plots
        bin_count, binedges, bin_nums = binned_statistic(
        curr_stat_df['Dist_X_Time'].to_numpy(), 
        curr_stat_df['d_ratio'].to_numpy(), bins = NUM_BINS, statistic= 'count')
        # print(min(bin_count))
        # print(bin_count)

        downsampled_df = []
        min_coverage = int(min(bin_count))
        num_unique_loci = []

        # #get the minimum number of distinct loci
        # for i in range(len(binedges)-1):
        #     bin_start = binedges[i]
        #     bin_end = binedges[i+1]
        #     curr_bin = curr_stat_df[(curr_stat_df['Dist_X_Time'] >= bin_start) & (curr_stat_df['Dist_X_Time'] <= bin_end)]
        #     unique_loci = set(list(curr_stat_df['Locus_1'].unique()) + list(curr_stat_df['Locus_2'].unique().tolist()))
        #     num_unique_loci.append(len(unique_loci))
        
        # min_unique_loci = int(min(num_unique_loci))

        # for i in range(len(binedges)-1):
        #     bin_start = binedges[i]
        #     bin_end = binedges[i+1]
        #     curr_bin = curr_stat_df[(curr_stat_df['Dist_X_Time'] >= bin_start) & (curr_stat_df['Dist_X_Time'] <= bin_end)]
        #     unique_loci = set(list(curr_stat_df['Locus_1'].unique()) + list(curr_stat_df['Locus_2'].unique().tolist()))
        #     sampled_loci = np.random.choice(list(unique_loci), min_unique_loci, replace = False)
        #     sample_df = curr_bin[curr_bin['Locus_1'].isin(sampled_loci) & (curr_bin['Locus_2'].isin(sampled_loci))]
        #     downsampled_df.append(sample_df)

        
        for i in range(len(binedges)-1):
            bin_start = binedges[i]
            bin_end = binedges[i+1]
            curr_bin = curr_stat_df[(curr_stat_df['Dist_X_Time'] >= bin_start) & (curr_stat_df['Dist_X_Time'] <= bin_end)]

            #sample segregating sites while the D' num is below the threshold
            sample_df = []
            sample_size = 0
            sampled_loci = set()
            unsampled_loci = set(list(curr_bin['Locus_1'].unique()) + list(curr_bin['Locus_2'].unique().tolist()))
            while sample_size < min_coverage and len(unsampled_loci) > 0:
                #sample another locus
                my_sample = np.random.choice(list(unsampled_loci))
                sampled_loci.add(my_sample)
                unsampled_loci.remove(my_sample)
                sample_df = curr_bin[curr_bin['Locus_1'].isin(sampled_loci) & (curr_bin['Locus_2'].isin(sampled_loci))]
                sample_size = len(sample_df)
                

            curr_bin = sample_df
            downsampled_df.append(curr_bin)

        downsampled_df = pd.concat(downsampled_df, ignore_index=True)
   
        sns.histplot(data = downsampled_df, stat = 'density', x = 'Dist_X_Time', color = 'blue', bins = NUM_BINS)

        plt.savefig(outDir + "fits/downsampled_hist"+ str(curr_float_rho) + ".png", dpi = 300)
        plt.close()


        #get the estimate and fit for the current dataset and sample size
        x_vals = downsampled_df['Dist_X_Time'].unique()
        coeffs, fit_dat = optimize.curve_fit(plne.neher_leitner, 
            downsampled_df['Dist_X_Time'], downsampled_df['d_ratio'],
            p0 = [0, 0.26, .0000439], maxfev = 10000)
        fit_vals = [plne.neher_leitner(x, coeffs[0], coeffs[1], coeffs[2])
                    for x in x_vals]

        #Bin the d' ratios so they are easier to view on the plots
        binned_rat, binedges, bin_nums = binned_statistic(
            downsampled_df['Dist_X_Time'].to_numpy(), 
            downsampled_df['d_ratio'].to_numpy(), bins = 100)

        sns.lineplot(y = binned_rat, x = binedges[1:], color = 'black')
        sns.lineplot(x = x_vals, y = fit_vals, color = 'red')
        plt.savefig(outDir + 'fits/fit_results' + str(curr_float_rho) + '.png')
        plt.close()

        estimate_df = pd.DataFrame([[coeffs[0], coeffs[1], coeffs[2], coeffs[1] * coeffs[2]]], columns = ['c0', 'c1', 'c2', 'Estimated_Rho'])

        estimate_df['Group'] = curr_iteration
        estimate_df['Sim_Rho'] = curr_rho
        all_estimate_df.append(estimate_df)


all_estimate_df = pd.concat(all_estimate_df, ignore_index=True)
print(all_estimate_df)
all_conf_ints = []

#Calculate the confidence intervals for the estimates with low + high VL
for name, group in all_estimate_df.groupby(['Group', 'Sim_Rho']):
    lower_conf = np.quantile(group['Estimated_Rho'], 0.025)
    mid_est = np.quantile(group['Estimated_Rho'], 0.5)
    upper_conf = np.quantile(group['Estimated_Rho'], 0.975)

    #append the confidence interval
    all_conf_ints.append([lower_conf, mid_est, upper_conf, name[0], name[1]])

all_conf_ints = pd.DataFrame(all_conf_ints, 
    columns=['lower_conf', 'est_rho', 'upper_conf', 'Group', 'Sim_Rho'])


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
for entry in all_conf_ints['Sim_Rho']:
    intRhoList.append(rhoDict[entry])
all_conf_ints['Sim_float_rho'] = intRhoList


########################## Plotting the accuracy of the estimates #############
#plot the estimates to show how accurate they are
sns.set(rc={'figure.figsize':(10,10)}, font_scale = 2, font = '')
sns.set_palette("tab10")
sns.set_style("white")

fig = sns.stripplot(x = 'Sim_Rho', y = 'est_rho', data = all_conf_ints, 
    jitter = True, color = 'k', s = 8, alpha = 0.3,
    order = [r"$2\times10^{-6}$", r"$10^{-5}$", r"$2\times10^{-5}$", r"$10^{-4}$", r"$2\times10^{-4}$" , r"$10^{-3}$"])
    #order = [r"$2\times10^{-6}$", r"$10^{-5}$", r"$2\times10^{-5}$", r"$10^{-4}$", r"$2\times10^{-4}$", r"$10^{-3}$"])

# distance across the "X" or "Y" stipplot column to span, in this case 40%
label_width = 0.4

#need to save the fig to force the axis to be drawn
plt.savefig(outDir + "selection_accuracy.png", dpi = 300)
            
for tick, text in zip(fig.get_xticks(), fig.get_xticklabels()):
    sample_name = text.get_text()  # "X" or "Y"

    #get the float value of rho corresponding with the tick
    rho_val = all_conf_ints[all_conf_ints['Sim_Rho'] == sample_name]
    rho_val = rho_val['Sim_float_rho'].unique()[0]

    # plot horizontal lines across the column, centered on the tick
    fig.plot([tick-label_width/2, tick+label_width/2], [rho_val, rho_val],
            lw=2, color='k')

plt.ylabel(r'Estimated Value of $\rho$')
plt.xlabel(r'Simulation Value of $\rho$')
plt.yscale('log')
plt.savefig(outDir + "downsampled_accuracy.png", dpi = 300)

#Need to convert the coefficient estimates into long form and then plot
# coefficients_for_plot = pd.melt(estimate_df, id_vars = ['Group', 'Sim_Rho'], value_vars = ['c0', 'c1', 'c2'])
# print(coefficients_for_plot)