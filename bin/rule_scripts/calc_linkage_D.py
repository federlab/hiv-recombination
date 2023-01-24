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
dataDir = dataDir.split('/')[:-2]
dataDir = "/".join(dataDir)

#The directories for our data
coCounts_dir = dataDir + "/numpy/"
loci_dir = dataDir + "/analysis/FilteredLoci"
genotype_dir = dataDir + "/analysis/FilteredGenotypes"

#Make the output directory if we need to 
if not os.path.exists(dataDir + "/linkage_D"):
    os.mkdir(dataDir + "/linkage_D")

#Load the segregating loci
segregatingLoci = pd.read_pickle(loci_dir)
genotype_df = pd.read_pickle(genotype_dir)

#One dataframe to put all of the linkage calcs in
all_resultsDF = []

#Loop through all of the timepoints files and get their information.
for curr_timepoint in segregatingLoci['timepoint'].unique():
    print(curr_timepoint, file = sys.stderr)

    #check the current segregating loci
    curr_seg = segregatingLoci[segregatingLoci['timepoint'] == curr_timepoint]
    curr_seg.reset_index(drop = True, inplace = True)

    #Get the current genotype information
    curr_genotypes = genotype_df[genotype_df['timepoint'] == curr_timepoint]
    curr_genotypes.reset_index(drop = True, inplace = True)

    #get the linkage calculation results
    stat_list, distList, supportList, resultsDF = r2.calculate_R2_df(
        curr_genotypes, curr_seg, statistic = 'D', saveData = True)
    
    resultsDF['timepoint'] = curr_timepoint
    all_resultsDF.append(resultsDF)

#Make the dataframe from all of our results
all_resultsDF = pd.concat(all_resultsDF, ignore_index=True)
all_resultsDF.to_pickle(dataDir + "/linkage_D/r2_and_D")
