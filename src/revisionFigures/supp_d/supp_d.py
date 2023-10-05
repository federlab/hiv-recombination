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
import slimUtil as slim_util
from scipy import optimize
from sklearn.metrics import mean_squared_error
from matplotlib import rcParams

DIST_TIME_MAX = 50000
NUM_BOOTSTRAPS = 1000
NUM_REPS = 200
NUM_GROUPS = 20
DI_THRESHOLD = 0.005

#This figure has the accuracy of the method with D values instead of D' values
dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/revision/supp_figs/supp_d/'
outDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/revision/supp_figs/supp_d/"

# dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/revision/supp_figs/supp_d/'
# outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/revision/supp_figs/supp_d/'

#Load the dataframes that will be plotted
all_conf_ints = pd.read_pickle(dataDir + "all_conf_ints_"+ \
                                        str(NUM_BOOTSTRAPS) + "_" + \
                                            str(NUM_GROUPS) + ".pkl")
all_stat_dfs = pd.read_pickle(dataDir + "all_stat_dfs_"+ \
                                        str(NUM_BOOTSTRAPS) + "_" + \
                                            str(NUM_GROUPS) + ".pkl")

all_conf_ints = slim_util.label_rho(all_conf_ints)
all_stat_dfs = slim_util.label_rho(all_stat_dfs)


########################## Plotting the accuracy of the estimates #############
#plot the estimates to show how accurate they are
rcParams.update({'figure.figsize':(3.5, 3), 'font.size': 8})


fig = sns.stripplot(x = 'Sim_Rho', y = 'est_rho', data = all_conf_ints,
    jitter = True, color = 'k', s = 3, alpha = 0.3,
    order = [r"$2\times10^{-6}$", r"$10^{-5}$", r"$2\times10^{-5}$", r"$10^{-4}$", r"$2\times10^{-4}$", r"$10^{-3}$"])
plt.savefig(outDir + "accuracy_D"  + str(NUM_BOOTSTRAPS) + ".png")
# distance across the "X" or "Y" stipplot column to span, in this case 40%
label_width = 0.4
                
for tick, text in zip(fig.get_xticks(), fig.get_xticklabels()):
    sample_name = text.get_text()  # "X" or "Y"
    print(sample_name)

    #get the float value of rho corresponding with the tick
    rho_val = all_conf_ints[all_conf_ints['Sim_Rho'] == sample_name]
    rho_val = rho_val['Sim_float_rho'].unique()[0]

    # plot horizontal lines across the column, centered on the tick
    fig.plot([tick-label_width/2, tick+label_width/2], [rho_val, rho_val],
            lw=1, color='k')

fig.set_ylabel(r'Estimated Recombination Rate ($\hat{\rho}$)')
fig.set_xlabel(r'Simulated Recombination Rate ($\rho$)')


plt.ylim(0.000001, 0.01)
plt.tight_layout()
plt.yscale('log')
plt.savefig(outDir + "accuracy_D" + str(NUM_BOOTSTRAPS) + ".png", dpi = 300)
plt.close()