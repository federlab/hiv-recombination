import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import os
import numpy as np
import pandas as pd
import estimation_util as est_util

DIST_TIME_MAX = 50000
RATIOS_PER_GROUP = 25000

NUM_BOOTSTRAPS = 1000
NUM_GROUPS = 8

#This file makes the estimated for the simulated data which are used to create
#the simulation figures for the manuscript

#For running on desktop
dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/'
outDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/paper/simulated_estimates_selection/"

# dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/slimDatasets/'
# outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/paper/simulated_estimates_selection/'
dataDirList = ['2022_09_07_MPL_1e-3/', '2022_09_07_MPL_1e-4/', '2022_09_07_MPL_1e-5/']
dataDirReps = [160, 3270, 5490]


#Get estimates for each of the selection rho values
for i in range(len(dataDirList)):
    curr_num_reps = dataDirReps[i]
    curr_dir = dataDir + dataDirList[i]
    curr_out = outDir + dataDirList[i]

    est_util.make_comparison_dataframes(curr_dir, curr_out, curr_num_reps, 
                                                NUM_GROUPS, NUM_BOOTSTRAPS, 
                                                DIST_TIME_MAX, RATIOS_PER_GROUP)