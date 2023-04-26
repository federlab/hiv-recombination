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

DIST_TIME_MAX = 50000
NUM_BOOTSTRAPS = 1000
NUM_GROUPS = 8

#Today I am going to prototype an analysis to determine how well we can
#discriminate between two different recombination rates that are a fixed
#distance apart.

dataDir =  "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/paper/simulated_estimates_selection/"
outDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/paper/supp_selection/"

# dataDir =  "/net/feder/vol1/home/evromero/2021_hiv-rec/results/paper/simulated_estimates_selection/"
# outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/paper/supp_selection/'
dataDirList = ['2022_09_07_MPL_1e-3/', '2022_09_07_MPL_1e-4/', '2022_09_07_MPL_1e-5/']

###############################################################################
#################################### Helper Functions #########################
###############################################################################

def label_rho(my_dataframe):
    #This function takes a dataframe and labels the rho values with floats
    #instead of strings so that they can be used as a continuous variable
    #To get the comparison, we need to convert the rho values to floats
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
    for entry in my_dataframe['Sim_Rho']:
        intRhoList.append(rhoDict[entry])
        newStringRho.append(rho_dict_fix_strings[entry])
    my_dataframe['Sim_float_rho'] = intRhoList
    my_dataframe['Sim_Rho'] = newStringRho

    return my_dataframe
###############################################################################
###############################################################################
###############################################################################
all_stat_dfs = []
all_conf_ints = []

#Get the D' ratios for all of the simulations
for curr_dir in dataDirList:
    curr_stat_df = pd.read_pickle(dataDir + curr_dir + "all_stat_dfs_"+ \
                                        str(NUM_BOOTSTRAPS) + "_" + \
                                            str(NUM_GROUPS) + ".pkl")
    all_stat_dfs.append(curr_stat_df)

    curr_conf_ints = pd.read_pickle(dataDir + curr_dir + "all_conf_ints_"+ \
                                        str(NUM_BOOTSTRAPS) + "_" + \
                                            str(NUM_GROUPS) + ".pkl")
    all_conf_ints.append(curr_conf_ints)

all_stat_dfs = pd.concat(all_stat_dfs, ignore_index=True)
all_stat_dfs = label_rho(all_stat_dfs)

all_conf_ints = pd.concat(all_conf_ints, ignore_index=True)
all_conf_ints = label_rho(all_conf_ints)

rho_dict_reverse = {r"$10^{-3}$" : 0.001,
                    r"$10^{-4}$" : 0.0001,
                    r"$2\times10^{-4}$" : 0.0002,
                    r"$10^{-5}$" : 0.00001,
                    r"$2\times10^{-5}$" : 0.00002,
                    r"$2\times10^{-6}$" : 0.000002}

########################## Plotting the accuracy of the estimates #############
#plot the estimates to show how accurate they are
rcParams.update({'figure.figsize':(4, 3), 'font.size': 6})


fig = sns.stripplot(x = 'Sim_Rho', y = 'est_rho', data = all_conf_ints, 
    jitter = True, color = 'k', s = 3, alpha = 0.5,
    order = [r"$2\times10^{-6}$", r"$10^{-5}$", r"$2\times10^{-5}$", r"$10^{-4}$", r"$2\times10^{-4}$", r"$10^{-3}$"])

# distance across the "X" or "Y" stipplot column to span, in this case 40%
label_width = 0.4

#need to save the fig to force the axis to be drawn
plt.savefig(outDir + "selection_accuracy"+ str(NUM_BOOTSTRAPS)+".png", dpi = 300)
            
for tick, text in zip(fig.get_xticks(), fig.get_xticklabels()):
    sample_name = text.get_text()  # "X" or "Y"

    rho_val = rho_dict_reverse[sample_name]

    if rho_val in [0.0002, 0.00002, 0.000002]:
        plt.annotate("N.D.", xy=(tick, rho_val), xytext= (tick - label_width/2, rho_val * 1.5))

    # plot horizontal lines across the column, centered on the tick
    fig.plot([tick-label_width/2, tick+label_width/2], [rho_val, rho_val],
            lw=1, color='k')

plt.ylabel(r'Estimated Recombination Rate ($\hat{\rho}$)')
plt.xlabel(r'Simulated Recombination Rate ($\rho$)')
plt.yscale('log')
plt.tight_layout()
plt.savefig(outDir + "selection_accuracy"+ str(NUM_BOOTSTRAPS)+".png", dpi = 300)
plt.close()

