import stat
import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
#for running on desktop
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import os
import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plot_neher as plne
import autocorrelation as autocorr
import zaniniUtil as zu
from scipy import optimize
from scipy import stats

# In this script I am repeating the viral load analysis, but I am resampling
# from the patients such that I leave one out in each sample
THRESHOLD = 0.2
DIST_TIME_MAX = 50000
GROUP_THRESHOLD_LIST = [10000, 50000, 100000, 200000]

#For running on desktop
dataDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/"
vlDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/"
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/paper/fig2/'

#Make the dataframe containg D' ratios
stat_df = zu.combine_drats(dataDir)
stat_df = zu.label_vl_drats(stat_df, vlDir)

#Filter the d' ratio dataframe
stat_df = stat_df[stat_df['Fragment'] != 'F5']
stat_df = stat_df[stat_df['Participant'] != 'p4']
stat_df = stat_df[stat_df['Time_Diff'].gt(50)]
stat_df = stat_df[stat_df['Monotonic'] == True]

################### Plot the viral loads ######################################
sns.ecdfplot(x = 'Ave_VL', data = stat_df)
plt.savefig(outDir + "vl_ecdf.jpg")
plt.close()

sns.histplot(x = 'Ave_VL', data = stat_df, cumulative = True)
plt.axvline(50000, color = 'black')
plt.ylabel("Total # of D' ratios")
plt.xscale('log')
plt.savefig(outDir + "vl_hist.jpg")
plt.close()