import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import os
import numpy as np
import pandas as pd
import estimation_util as est_util

DIST_TIME_MAX = 50000
LOCI_PER_GROUP = 1000
THRESHOLD = 0.005

NUM_BOOTSTRAPS = 1000
NUM_REPS = 200
NUM_GROUPS = 20

#This file makes the estimated for the simulated data which are used to create
#the simulation figures for the manuscript

#For running on desktop
dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/paper/simulated_estimates_d/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/revision/supp_figs/supp_d/'

dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/paper/simulated_estimates_d/'
outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/revision/supp_figs/supp_d/'

est_util.make_comparison_dataframes(dataDir, outDir, NUM_REPS, 
                                            NUM_GROUPS, NUM_BOOTSTRAPS, 
                                            DIST_TIME_MAX, THRESHOLD = THRESHOLD,
                                            LOCI_PER_GROUP = LOCI_PER_GROUP,
                                            PREMADE_DF = dataDir + "all_stat_dfs_" + \
                               str(NUM_BOOTSTRAPS) + "_" + str(NUM_GROUPS) \
                                + ".pkl") 