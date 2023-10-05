import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import estimation_util as est_util
import slimUtil as slim_util
from matplotlib import rcParams
from matplotlib.lines import Line2D
from scipy.stats import binned_statistic


#This figure is in response to Editor comment 1.1
#First, as you mention in lines 102-106, linkage disequilibrium results from 
#new mutations that occur in particular genetic background and then 
#disassociate with recombination. Then, what you may expect from a within-host
#population of HIV is an equilibrium level of LD due to recurrent mutation and
#recombination, particularly when one anticipates a characteristic genealogy
#(phylogeny) of within-host HIV (for example, Figure 1I of Grenfell et al. 
#2004) that shows continuous turn-over of viral lineages. I therefore expect
#that what is driving the relationship between D' ratio and distance*time in
#Figure 2 is mostly distance (at a given time). I wonder whether you actually
#see D' decrease over time while distance is fixed.


NUM_BOOTSTRAPS = 1000
NUM_REPS = 500
NUM_GROUPS = 20

#The rho values to show in the accuracy panel
ACC_LIST = [r"$2\times10^{-6}$",
                        r"$5\times10^{-6}$",
                        r"$10^{-5}$",
                        r"$2\times10^{-5}$",
                        r"$5\times10^{-5}$",
                        r"$10^{-4}$",
                        r"$2\times10^{-4}$",
                        r"$5\times10^{-4}$",
                        r"$10^{-3}$"]

dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/revision/simulated_estimates_neutral/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/revision/reviewer_dist/'

######################### Configure Plot Settings #############################
#plot the estimates to show how accurate they are
params = {'figure.figsize':(7.5, 3), 'axes.labelsize': 6,'axes.titlesize':6,  
          'legend.fontsize': 6, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
          'legend.title_fontsize': 6}
linewidth = 1
markerSize = 3
rcParams.update(params)


###############################################################################
###############################################################################
###############################################################################
#Load the dataframes that will be plotted
all_stat_dfs = pd.read_pickle(dataDir + "all_stat_dfs_"+ \
                                        str(NUM_BOOTSTRAPS) + "_" + \
                                            str(NUM_GROUPS) + ".pkl")

all_stat_dfs = slim_util.label_rho(all_stat_dfs)
all_stat_dfs['Dist'] = all_stat_dfs['Locus_2'] - all_stat_dfs['Locus_1']
print(all_stat_dfs.columns)

sns.histplot(data=all_stat_dfs, x='Dist', hue='Time_Diff', alpha = 0.3)
plt.savefig(outDir + 'Dist_hist.png')