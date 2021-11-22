import sys
#for cluster run
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import neher
import zaniniUtil as zu

#directories for cluster run
dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/analysis/'
viralLoadDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/'
outData = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/analysis/neher_downsampled/'

#Today I am going to conduct the neher and leitner analysis
#But I am going to downsample multinomially at each timepoint
#I need to do this, save the data, then match up the timepoints with viral loads.
#Note: we are not filtering the haplotypes for this

#Next we need to set some things for the run
fragment_list = ['F1','F2', 'F3', 'F4', 'F5','F6']
par_list = ['p1', 'p2', 'p3', 'p5', 'p6', 'p8', 'p9', 'p10', 'p11']
par_list =['p1']
fragment_list = ['F1']
available_files_hap = os.listdir(dataDir + "haplotype_dfs/")
available_files_seg = os.listdir(dataDir + "segregating_Loci/")


NUM_VL_GROUPS = 5
BINWIDTH = 100
MIN_BIN = 0
MAX_BIN = 700
CUTOFF = 0.03
#The frequency a haplotype must have to trigger a test success
SUCCESSFILTER = 0
RUNNAME = str(CUTOFF) +  "filtering_success_"  + str(SUCCESSFILTER) +"_binnedVL.png"

#We can start by getting the viral load data
viralLoadData = zu.make_viral_load_df(viralLoadDir)
#We can also load our dataframe for converting days to timepoints
dayToTime = pd.read_csv(dataDir + "DaysToTimepoints.csv")
labeledLoads = []

#We need to add a timepoint column to our viral load data
for curr_par in par_list:
    curr_data = viralLoadData[viralLoadData['Participant'] == curr_par]
    curr_days = dayToTime[dayToTime['Participant'] == curr_par]
    labels = []
    #get the timepoint for the day
    for index, row in curr_data.iterrows():
        date = int(row['Days from infection'])
        currLabel = curr_days[curr_days[' Day'] == date]
        # print("****************", file = sys.stderr)
        # print(curr_par, file = sys.stderr)
        # print(currLabel, file = sys.stderr)
        # print(date, file = sys.stderr)
        labels.append(currLabel[' Timepoint'].tolist()[0])
    labeled_current = curr_data.copy()
    labeled_current['Timepoint'] = labels
    labeledLoads.append(labeled_current)
labeledLoads = pd.concat(labeledLoads)
viralLoadData = labeledLoads

#make a dataframe to store our data to plot
all_frequencies_patients = []

#loop through the patients and get their results from the neher analysis
for curr_par in par_list:
    for curr_fragment in fragment_list:
        par_frag = curr_par +  "_" + curr_fragment
        haplotype_file = "haplotypes_" + par_frag + ".pkl"
        loci_file = "segregating_" + par_frag + ".pkl"
        
        #check if there are  files for this
        if haplotype_file not in available_files_hap or loci_file not in available_files_seg:
            continue

        #First we need to load in the dataframes for each patient + fragment.
        haplotype_df = pd.read_pickle(dataDir + "haplotype_dfs/" + haplotype_file)
        segregating_Loci = pd.read_pickle(dataDir + "segregating_Loci/" + loci_file)

        #now we can perform neher analysis on this dataframe
        #first we will filter out genotypes with alleles at less than 3% frequency
        haplotype_df = zu.filter_genotype_df(haplotype_df, segregating_Loci, allele_cutoff = CUTOFF, hap_cutoff = SUCCESSFILTER)
        if haplotype_df.empty:
            continue

        print(haplotype_df, file = sys.stderr)
        quit()
