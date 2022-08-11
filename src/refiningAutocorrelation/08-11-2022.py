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

THRESHOLD = 0.2

dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_08_10_MPL/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_08_10_MPL/'

estimate_df = [] 
seg_hist_data = []
all_stat_df = []

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

    #make a list of all the segregating loci
    all_seg_loc_1 = set(stat_df["Locus_1"].unique())
    all_seg_loc_2 = set(stat_df["Locus_2"].unique())
    all_seg_loc = all_seg_loc_1.union(all_seg_loc_2)
    all_seg_loc = np.array(list(all_seg_loc))
    seg_hist_data.append([len(all_seg_loc), sim_rho])

    stat_df['sim_rho'] = sim_rho
    stat_df['rep'] = rep
    all_stat_df.append(stat_df)

seg_hist_data = pd.DataFrame(seg_hist_data, columns = ['num_seg_loci', 'sim_rho'])
seg_hist_data = seg_hist_data.groupby('sim_rho').num_seg_loci.agg('sum')
all_stat_df = pd.concat(all_stat_df, ignore_index=True)
print(seg_hist_data)
sim_rho_list = all_stat_df['sim_rho'].unique()

for curr_rho in sim_rho_list:
    stat_df = all_stat_df[all_stat_df['sim_rho'] == curr_rho]

    #get the estimate and fit for the current dataset and sample size
    x_vals = stat_df['Dist_X_Time'].unique()
    coeffs, fit_dat = optimize.curve_fit(plne.neher_leitner, 
        stat_df['Dist_X_Time'].to_list(), stat_df['d_ratio'].to_list(),
        p0 = [0, 0.26, .0000439], maxfev = 10000)
    fit_vals = [plne.neher_leitner(x, coeffs[0], coeffs[1], coeffs[2])
                for x in x_vals]

    #Bin the d' ratios so they are easier to view on the plots
    binned_rat, binedges, bin_nums = binned_statistic(
        stat_df['Dist_X_Time'].to_list(), 
        stat_df['d_ratio'].to_list(), bins = 100)

    estimate_df.append([coeffs[0], coeffs[1], coeffs[2], 
        coeffs[1] * coeffs[2], curr_data, curr_rho, 
        curr_data])  

    sns.lineplot(x = binedges[:-1], y = binned_rat, color = 'tab:blue')
    sns.lineplot(x = x_vals, y = fit_vals, color = 'tab:orange')
    plt.xlabel(r'Distance X Time (bp/generation)')
    plt.ylabel("D\' Ratio")


    plt.tight_layout()
    plt.savefig(outDir + curr_rho + 'fit.jpg')
    plt.close()

estimate_df = pd.DataFrame(estimate_df, columns=["C0", "C1", "C2",
                     "Est_Rho", 'Dataset', 'Sim_Rho', 'data'] )

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

fig = sns.stripplot(x = 'Sim_Rho', y = 'Est_Rho', data = estimate_df, 
    jitter = True, color = 'k', s = 8,
    order = [r"$2\times10^{-6}$", r"$10^{-5}$", r"$2\times10^{-5}$", r"$10^{-4}$", r"$2\times10^{-4}$", r"$10^{-3}$"])

# distance across the "X" or "Y" stipplot column to span, in this case 40%
label_width = 0.4

            
for tick, text in zip(fig.get_xticks(), fig.get_xticklabels()):
    sample_name = text.get_text()  # "X" or "Y"

    #get the float value of rho corresponding with the tick
    rho_val = estimate_df[estimate_df['Sim_Rho'] == sample_name]
    rho_val = rho_val['Sim_int_rho'].unique()[0]

    # plot horizontal lines across the column, centered on the tick
    fig.plot([tick-label_width/2, tick+label_width/2], [rho_val, rho_val],
            lw=2, color='k')


plt.xlabel(r'Simulation Value of $\rho$')
plt.ylabel(r'Estimated Value of $\rho$')
plt.ylim(0.000001, 0.01)
plt.yscale('log')
plt.tight_layout()

plt.savefig(outDir + "neut_test.jpg")
plt.close()

sns.histplot(seg_hist_data, bins = 20)
plt.savefig(outDir + "neut_hist.jpg")
plt.close()

