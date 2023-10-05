import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import pandas as pd
import slimUtil as slim_util

#I already compressed the selection data so we will be using some of the vNE data which has the exact 
#same population parameters as the selection data (the 10^4 subset of it)
dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2023_08_28_selection_vNE/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/revision/simulated_linkage_sel/'

dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2023_08_28_selection_vNE/'
outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/revision/simulated_linkage_sel/'

AF_THRESH = (0.2, 0.8)
AF_THRESH_STR = ('20', '80')
DIST_THRESH = 400

all_D_vals_df = slim_util.make_dprime_df(dataDir, 
                    allele_freq_threshes = AF_THRESH, distance_thresh = DIST_THRESH,
                    four_hap_thresh = True, migration = False, rho_list = ['1e-04', '1e-05'],
                    rep_limit = 250)

fileTags = AF_THRESH_STR[0] + '_' + AF_THRESH_STR[1] + '_' + str(DIST_THRESH) + '.pkl'
all_D_vals_df.to_pickle(outDir + 'all_D_vals_df' + fileTags)