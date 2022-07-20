import sys
#for running on cluster
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
import os
import neher
import numpy as np
import pandas as pd
import zaniniUtil as zu

#This script takes in the cocount arrays and makes a dataframe of 
#filtered genotypes
#I am going to use the same filtering as our previous analysis
#These parameters are used when processing the data
CUTOFF = 0.03
SUCCESS = 0.01
RUNNAME = str(CUTOFF) + '_' + str(SUCCESS) +  ""

dataDir = snakemake.output[0]
print(dataDir, file = sys.stderr)
dataDir = dataDir.split('/')[:-2]
dataDir = "/".join(dataDir)

timepoint_dir = dataDir + "/numpy/"


#make a dictionary for the timepoint labels
timepoint_df = pd.read_csv(dataDir + '/timepoint_info.tsv', sep = ' ',
                header= None, names = ['name', 'generation'], index_col = False)

#a dataframe to save all of the timepoints in
all_timepoint_genotypes = []

#Loop through all of the files and get their information.
for currfile in os.listdir(timepoint_dir):

    #make sure the file has the right type
    if currfile.split('.')[-1] != "npy":
        continue

    #get the timepoint label
    timepoint = currfile.split('_')[-1]
    timepoint = timepoint.split('.')[0]

    #load the array
    coCounts_arr = np.load(timepoint_dir  + currfile)
    frag_len = coCounts_arr.shape[-1]

    #find the segregating sites
    segregatingLoci = zu.find_segregating_diagonal(coCounts_arr, all_seg = True)  
    segregatingLoci['timepoint'] = timepoint_df[timepoint_df['name'] == int(timepoint[-1])]['generation'].tolist()[0]
    segregatingLoci['frag_len'] = frag_len

    #make our dataframe of genotypes
    genotype_df = zu.make_genotype_df(segregatingLoci, coCounts_arr)

    #filter the dataframe 
    genotype_df['timepoint'] = timepoint_df[timepoint_df['name'] == int(timepoint[-1])]['generation'].tolist()[0]
    filtered_genotypes = zu.filter_genotype_df(genotype_df, segregatingLoci, allele_cutoff= CUTOFF, hap_cutoff= SUCCESS)
    all_timepoint_genotypes.append(filtered_genotypes)


all_timepoint_genotypes = pd.concat(all_timepoint_genotypes)
all_timepoint_genotypes.to_pickle(dataDir + '/analysis/FilteredGenotypes')
print(all_timepoint_genotypes, file = sys.stderr)