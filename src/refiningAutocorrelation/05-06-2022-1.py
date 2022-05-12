import os
import pandas as pd
#I am going to use this file make the directory structure necessary for 
#running the snakemake pipeline on the Zanini data.

dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini/'
outDataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/'

#first make a directory for each participant and fragment
for curr_file in os.listdir(dataDir + 'snpPairs'):
    #ignore DS_Store etc
    if curr_file.split('.')[-1] != 'npy':
        continue

    curr_file = curr_file.split('.')[0]
    par = curr_file.split('_')[1]
    frag = curr_file.split('_')[4]
    par_frag  = par + "_" + frag
    
    if not os.path.exists(outDataDir + par_frag):
        os.mkdir(outDataDir + par_frag)

#get the master timepoint file
timepointFile = dataDir + "analysis/DaysToTimepoints.csv"

timepoint_df = pd.read_csv(timepointFile)
timepoint_df['timepoint'] = 0.5 * timepoint_df[' Day']


#write a timepoint file into the folder for each participant and fragment
for curr_dir in os.listdir(outDataDir):
    #skip the DS store and other hidden files
    if curr_dir[0] == '.':
        continue
    par = curr_dir.split("_")[0]
    frag = curr_dir.split("_")[1]
    curr_info = timepoint_df[timepoint_df['Participant'] == par]
    curr_info = pd.concat([curr_info[' Timepoint'], curr_info['timepoint']], axis = 1)
    curr_info.set_index(' Timepoint', inplace= True)
    curr_info.index.name = None
    curr_info.to_csv( outDataDir + curr_dir +'/timepoint_info.tsv', sep = ' ', header = False)
