import sys
#for running on cluster
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
import os
import numpy as np
import pandas as pd
import autocorrelation as autocorr

THRESHOLD = 0.2

#This script uses the coCounts_arr and list of segregatingLoci
#It takes them and saves a dataframe with the R^2 and D statistics
dataDir = snakemake.output[0]
print(dataDir, file = sys.stderr)
dataDir = dataDir.split('/')[:-2]
dataDir = "/".join(dataDir)

linkage_file = dataDir + "/linkage/r2_and_D"
d_ratio_out =  dataDir + "/linkage/d_ratio"

stat_df = autocorr.calculate_d_ratios(linkage_file, THRESHOLD)
stat_df.to_pickle(d_ratio_out)
print(os.listdir(dataDir + "/linkage/"))
if 'd_ratio' not in os.listdir(dataDir + "/linkage/"):
    print(stat_df, file = sys.stderr)
    stat_df.to_pickle(d_ratio_out)