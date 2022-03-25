import sys
from tabnanny import verbose
# sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
#for running on desktop
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import zaniniUtil as zu
import r2Analysis as r2

test_loci = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_02_24/mu1e-05_rho1e-04_Ne10000_M800_rep1/analysis/FilteredLoci"
test_loci = pd.read_pickle(test_loci)

coCounts_arr = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_02_24/mu1e-05_rho1e-04_Ne10000_M800_rep1/numpy/slim_formatted_t1.npy"
coCounts_arr = np.load(coCounts_arr)

timepoint_df = pd.read_csv('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_02_24/mu1e-05_rho1e-04_Ne10000_M800_rep1/timepoint_info.tsv', sep = ' ',
                header= None, names = ['name', 'generation'], index_col = False)

timepoint = 1
generation = timepoint_df[timepoint_df['name'] == timepoint]
generation = generation['generation'].tolist()[0]

#we can use our functions to save the data
stat_list, distList, supportList, resultsDF = r2.calculate_R2_pairCounts(coCounts_arr, test_loci, statistic = 'D', saveData = True)
print(resultsDF)
