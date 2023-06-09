import sys
#for running on cluster
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
# #for running on cluster
# sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import os
import numpy as np
import pandas as pd
import autocorrelation as autocorr

THRESHOLD = -float('inf')

#This script uses the linkage data and outputs calculated D' ratios
dataDir = snakemake.input[0]
print(dataDir, file = sys.stderr)
dataDir = dataDir.split('/')[:-2]
dataDir = "/".join(dataDir)

outDataDir = snakemake.output[0]
outDataDir = outDataDir.split('/')[:-2]
outDataDir = "/".join(outDataDir)

linkage_file = dataDir + "/linkage_D/r2_and_D"
d_ratio_out =  outDataDir + "/linkage_D/d_ratio_three_haps"
r2_and_D = pd.read_pickle(linkage_file)


stat_df = autocorr.calculate_d_ratios(linkage_file, threshold = THRESHOLD, four_haps = False, stat = 'd_val')
stat_df.to_pickle(d_ratio_out)
