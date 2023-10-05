import sys
import pandas as pd
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import seaborn as sns
import numpy as np
import estimation_util as est_util
import slimUtil as slim_util
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

DIST_TIME_MAX = 50000
NUM_BOOTSTRAPS = 1000
NUM_REPS = 5000
NUM_GROUPS = 5

SAMPLE_RHO_LIST = [ r"$10^{-5}$", r"$10^{-4}$"]

#In this file I am going to be investigating the effects of decreased sampling depth
dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2023_09_12_selection_samp/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/revision/simulated_estimates_sampling/'

dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2023_09_12_selection_samp/'
outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/revision/simulated_estimates_sampling/'

est_util.samp_make_comparison_dataframes(dataDir, outDir, NUM_REPS, 
                                            NUM_GROUPS, NUM_BOOTSTRAPS, 
                                            DIST_TIME_MAX, DOWNSAMPLE = False) 