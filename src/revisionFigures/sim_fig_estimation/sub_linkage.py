import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import pandas as pd
import slimUtil as slim_util

dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2023_09_01_selection_sub/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/revision/simulated_linkage_sub/'

dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2023_09_01_selection_sub/'
outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/revision/simulated_linkage_sub/'

AF_THRESH = (0.2, 0.8)
AF_THRESH_STR = ('20', '80')
DIST_THRESH = 400

all_D_vals_df = slim_util.make_dprime_df(dataDir, 
                    allele_freq_threshes = AF_THRESH, distance_thresh = DIST_THRESH,
                    four_hap_thresh = True, migration = True)

fileTags = AF_THRESH_STR[0] + '_' + AF_THRESH_STR[1] + '_' + str(DIST_THRESH) + '.pkl'
all_D_vals_df.to_pickle(outDir + 'all_D_vals_df' + fileTags)