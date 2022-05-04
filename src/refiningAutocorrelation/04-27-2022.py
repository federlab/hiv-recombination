import sys
# sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
#for running on desktop
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plot_neher as plne
import autocorrelation as autocorr
from scipy import optimize
from scipy import stats

#This script just makes a basic plot of autoccorelation estimates for simulated data against
#the actual values in the simulation
#from this analysis, it looks like fitting on the sampled data doesn't quite work
THRESHOLD = 0.2
DIST_TIME_MAX = 50000
NUM_BINS = 25
NUM_PER_BIN = 500 #number of observations per bin (sampled with replacement)

#make a set of bins and sample similar amounts in every bin
bin_starts = range(0, DIST_TIME_MAX, int(DIST_TIME_MAX/NUM_BINS))
bin_edges = [(x, x + int(DIST_TIME_MAX/NUM_BINS)) for x in bin_starts]

#For running on cluster
dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_04_20/'
outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_04_20/'

dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_04_20/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_04_20/04-26-2022/'

estimate_df = [] 

for curr_data in os.listdir(dataDir):
    print(curr_data)

    #only get the data directories, not hidden files
    if curr_data[0] == '.':
        continue
    run_info = curr_data.split('_')
    sim_rho = run_info[1]
    sim_rho = sim_rho[3:]
    rep = run_info[-1]

    #make a place to store our output
    currOut = outDir + curr_data
    if not os.path.exists(currOut):
        os.mkdir(currOut)

    linkage_file = dataDir + curr_data + "/linkage/r2_and_D"
    d_ratio_file = dataDir + curr_data + "/linkage/d_ratio"

    if os.path.exists(d_ratio_file):
        stat_df = pd.read_pickle(d_ratio_file)
    else:
        stat_df = autocorr.calculate_d_ratios(linkage_file, THRESHOLD)
        stat_df.to_pickle(d_ratio_file)

    stat_df['dist'] = stat_df['Dist_X_Time']/stat_df['Time_Diff']

    bin_obs = []
    sampled_data = []
    data_aves = []
    bin_edges = bin_edges[2:]
    for curr_bin in bin_edges:
        stat_df_curr = stat_df[stat_df['Dist_X_Time'].between(curr_bin[0], curr_bin[1])]

        x_vals = stat_df_curr['Dist_X_Time'].unique()
        bin_obs.append(len(stat_df_curr))
        sampled_data.append(stat_df_curr.sample(n=500, replace = True))

        data_aves.append([curr_bin, np.mean([curr_bin[0], curr_bin[1]]) ,np.mean(stat_df_curr['d_ratio']), np.std(stat_df_curr['d_ratio'])/np.sqrt(len(stat_df_curr))])
        
    sampled_data = pd.concat(sampled_data, ignore_index=True)
    data_aves = pd.DataFrame(data_aves, columns = ['bin', 'Dist_X_Time', 'd_ratio', 'error'])
    # print(len(stat_df_curr))
    #Fit the data and add our estimate to the dataframe
    coeffs, fit_dat = optimize.curve_fit(plne.neher_leitner, sampled_data['Dist_X_Time'], sampled_data['d_ratio'], p0 = [0, 0.26, .0000439])
    fit_vals = [plne.neher_leitner(x, coeffs[0], coeffs[1], coeffs[2]) for x in x_vals]
    # sns.scatterplot(x = 'Dist_X_Time', y = 'd_ratio', data = sampled_data, alpha = 0.05, hue = 'd_i')
    sns.scatterplot(x = 'Dist_X_Time', y = 'd_ratio', data = data_aves)
    plt.errorbar(x = data_aves['Dist_X_Time'], y = data_aves['d_ratio'], yerr = data_aves['error'], xerr = None, ls = 'none', ecolor = 'gray')
    # sns.lineplot(x = 'Dist_X_Time', y = 'd_ratio', data = sampled_data, estimator = np.mean)
    plt.savefig(currOut + "/auto_plot_binned" + str(DIST_TIME_MAX) +".jpg")
    plt.close()


    #statistic, bin_edges, bin_number = stats.binned_statistic()

    # #what if we test fitting things just to the first 5k and then getting the saturation from a line at the later timepoints
    # stat_df1k = stat_df[stat_df['Dist_X_Time'].between(0, 1000)]
    # coeffs1k, fit_dat1k = optimize.curve_fit(plne.neher_leitner, stat_df1k['Dist_X_Time'], stat_df1k['d_ratio'], p0 = [0, 0.26, .0000439])
    # stat_df20_50k = stat_df[stat_df['Dist_X_Time'].between(20000, 50000)]
    # coeffs20_50k, fit_dat1k = optimize.curve_fit(plne.neher_leitner, stat_df1k['Dist_X_Time'], stat_df1k['d_ratio'], p0 = [0, 0.26, .0000439])
    # print(coeffs1k[2] * coeffs20_50k[1])


    estimate_df.append([coeffs[0], coeffs[1], coeffs[2], coeffs[1] * coeffs[2], curr_data, sim_rho, DIST_TIME_MAX, 'curve', curr_data])  
    # print(estimate_df)
    # quit()


estimate_df = pd.DataFrame(estimate_df, columns=["C0", "C1", "C2", "Est_Rho", 'Dataset', 'Sim_Rho', 'distTimeMax', 'style', 'data'] )


############################# Plotting Estimate Accuracy ######################
# #Plot our estimates against each other 
#make the rho values ints rather than strings
rhoDict = {"0.001" : 0.001,
            "1e-04" : 0.0001,
            "2e-04" : 0.0002,
            "1e-05" : 0.00001,
            "2e-05" : 0.00002,
            "2e-06" : 0.000002}

intRhoList = []
for entry in estimate_df['Sim_Rho']:
    intRhoList.append(rhoDict[entry])
estimate_df['Sim_int_rho'] = intRhoList

ax = sns.stripplot(x = 'Sim_Rho', y = 'Est_Rho', data = estimate_df, jitter = True,
    order = ["0.001", "2e-04", "1e-04", "2e-05", "1e-05", "2e-06"])

# distance across the "X" or "Y" stipplot column to span, in this case 40%
label_width = 0.4

for tick, text in zip(ax.get_xticks(), ax.get_xticklabels()):
    sample_name = text.get_text()  # "X" or "Y"

    # calculate the median value for all replicates of either X or Y
    rho_val = rhoDict[sample_name]

    # plot horizontal lines across the column, centered on the tick
    ax.plot([tick-label_width/2, tick+label_width/2], [rho_val, rho_val],
            lw=2, color='k')

plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
plt.xlabel("Simulation Value of Rho")
plt.ylabel("Estimated Value of Rho")
plt.ylim(0.000000001, 0.1)
plt.yscale('log')
plt.tight_layout()
plt.savefig(outDir + "autocorr_comparedEstimates_stripplot_" + str(DIST_TIME_MAX) + "_" + str(NUM_BINS) + "_" + str(NUM_PER_BIN)+".jpg")
plt.close()