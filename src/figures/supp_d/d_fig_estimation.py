import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import os
import numpy as np
import pandas as pd
import estimation_util as est_util

DIST_TIME_MAX = 50000
RATIOS_PER_GROUP = 25000
THRESHOLD = 0.005

NUM_BOOTSTRAPS = 1000
NUM_REPS = 200
NUM_GROUPS = 20

#This file makes the estimated for the simulated data which are used to create
#the simulation figures for the manuscript

#For running on desktop
dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_10_03_neutral/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/paper/simulated_estimates_d/'

# dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_10_03_neutral/'
# outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/paper/simulated_estimates_d/'

est_util.make_comparison_dataframes(dataDir, outDir, NUM_REPS, 
                                            NUM_GROUPS, NUM_BOOTSTRAPS, 
                                            DIST_TIME_MAX, RATIOS_PER_GROUP,
                                            THRESHOLD, True) 