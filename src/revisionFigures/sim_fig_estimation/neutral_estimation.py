import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import estimation_util as est_util


DIST_TIME_MAX = 50000
LOCI_PER_GROUP = 1000
NUM_BOOTSTRAPS = 1000
NUM_REPS = 100
NUM_GROUPS = 20


#This script runs the bootstrapped estimation for the manuscript
dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/revision/simulated_estimates_neutral/'
outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/revision/simulated_estimates_neutral/'


est_util.make_comparison_dataframes(dataDir, outDir, NUM_REPS, 
                                            NUM_GROUPS, NUM_BOOTSTRAPS, 
                                            DIST_TIME_MAX, LOCI_PER_GROUP = LOCI_PER_GROUP, 
                                            PREMADE_DF = dataDir + "all_stat_dfs_" + \
                               str(NUM_BOOTSTRAPS) + "_" + str(NUM_GROUPS) \
                                + ".pkl") 