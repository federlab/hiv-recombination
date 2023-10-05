import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import estimation_util as est_util


DIST_TIME_MAX = 50000
NUM_BOOTSTRAPS = 1000
NUM_REPS = 250
NUM_GROUPS = 5


#This script runs the bootstrapped estimation for the manuscript
dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2023_09_01_selection_sub/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/revision/simulated_estimates_sub/'

dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2023_09_01_selection_sub/'
outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/revision/simulated_estimates_sub/'

est_util.mig_make_comparison_dataframes(dataDir, outDir, NUM_REPS, 
                                            NUM_GROUPS, NUM_BOOTSTRAPS, DIST_TIME_MAX,
                                            RATIOS_PER_GROUP = None, DOWNSAMPLE = False) 

