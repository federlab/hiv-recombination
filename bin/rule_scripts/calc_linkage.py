import sys
#for running on cluster
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
import os
import numpy as np
import pandas as pd
import r2Analysis as r2

#This script uses the coCounts_arr and list of segregatingLoci
#It takes them and saves a dataframe with the R^2 and D statistics
dataDir = snakemake.output[0]
print(dataDir, file = sys.stderr)
dataDir = dataDir.split('/')[:-2]
dataDir = "/".join(dataDir)


#make a dictionary for the timepoint labels
timepoint_df = pd.read_csv(dataDir + '/timepoint_info.tsv', sep = ' ',
                header= None, names = ['name', 'generation'], index_col = False)
print(timepoint_df, file = sys.stderr)
#The directories for our data
coCounts_dir = dataDir + "/numpy/"
loci_dir = dataDir + "/analysis/FilteredLoci"

#Load the segregating loci
segregatingLoci = pd.read_pickle(loci_dir)

#One dataframe to put all of the linkage calcs in
all_resultsDF = []

#Loop through all of the timepoints files and get their information.
for currfile in os.listdir(coCounts_dir):
    #make sure the file has the right type
    if currfile.split('.')[-1] != "npy":
        continue

    #get the timepoint label
    timepoint = currfile.split('_')[-1]
    timepoint = timepoint.split('.')[0]
    timepoint = int(timepoint[1:])
    generation = timepoint_df[timepoint_df['name'] == timepoint]
    generation = generation['generation'].tolist()[0]

    #check the current segregating loci
    curr_seg = segregatingLoci[segregatingLoci['timepoint'] == generation]
    curr_seg.reset_index(drop = True, inplace = True)

    #load the array
    coCounts_arr = np.load(coCounts_dir  + currfile)

    #get the linkage calculation results
    stat_list, distList, supportList, resultsDF = r2.calculate_R2_pairCounts(
        coCounts_arr, curr_seg, statistic = 'D', saveData = True)
    
    resultsDF['timepoint'] = generation
    all_resultsDF.append(resultsDF)

#Make the dataframe from all of our results
all_resultsDF = pd.concat(all_resultsDF, ignore_index=True)
all_resultsDF.to_pickle(dataDir + "/linkage/r2_and_D")


