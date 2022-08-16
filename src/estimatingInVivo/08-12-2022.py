import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
#for running on desktop
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plot_neher as plne
import zaniniUtil as zu
from scipy import stats
from matplotlib import rcParams

#Today I am going to try to estimate the number of segregating loci in the
#in vivo data. Specifically I want to know how many segregating loci and pairs
#of sites go into the estimation for each group

GROUP_THRESHOLD_LIST = [7500, 10000, 25000, 50000, 100000, 200000]
NUM_BOOTSTRAPS = 1000

#For running on desktop
dataDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/"
vlDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/"
outDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/zanini/08-12-2022/"

#Make the dataframe containg D' ratios
stat_df = zu.combine_drats(dataDir)
stat_df = zu.label_vl_drats(stat_df, vlDir)

stat_df = stat_df[stat_df['Fragment'] != 'F5']
stat_df = stat_df[stat_df['Time_Diff'].gt(50)]

#Only get values for which the viral load doesn't leave the boundaries of the 
#initial and final measurments
stat_df = stat_df[stat_df['Monotonic'] == True]

#make a place to store the segregating loci counts
seg_loci_df = []
paired_loci_df = []

for curr_thresh in GROUP_THRESHOLD_LIST:
    #Get the dataframe for everyone except the current participant
    stat_df['High_VL'] = stat_df['Ave_VL'].gt(curr_thresh)
    high_vl_seg = 0
    low_vl_seg = 0

    #count the number of segregating loci in each group
    for name,group in stat_df.groupby(by = ['Participant', 'Locus_1', 'Locus_2', 'High_VL']):
        if group['High_VL'].unique()[0] == True:
            high_vl_seg += 1
        else:
            low_vl_seg += 1
    
    #add our current counts to the dataframe
    seg_loci_df.append([curr_thresh, high_vl_seg, low_vl_seg])

    #count the number of pairs of segregating loci in each group
    paired_loci_df.append([curr_thresh, len(stat_df[stat_df['High_VL'] == True]), len(stat_df[stat_df['High_VL'] == False])])

seg_loci_df = pd.DataFrame(seg_loci_df, columns = ['Threshold', 'High_VL', 'Low_VL'])
paired_loci_df = pd.DataFrame(paired_loci_df, columns = ['Threshold', 'High_VL', 'Low_VL'])
print('Counts of pairs')
print(paired_loci_df)
print('Counts of segregating loci')
print(seg_loci_df)
