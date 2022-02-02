import sys
#for running on cluster
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
# #for running on desktop
# sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import neher
import zaniniUtil as zu

#This file takes in a simulated dataset from slim
#It uses processes and filters the segregating sites then saves the filtered
#I am going to use the same filtering as our previous analysis
#These parameters are used when processing the data
CUTOFF = 0.03
SUCCESS = 0.01
RUNNAME = str(CUTOFF) + '_' + str(SUCCESS) +  ""

# #For running on desktop
# dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/'

#For running on cluster
dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/slimDatasets/'

enclosing_dir = '2022_01_25/'


#Today I am going to be using the data that Alison simulated using SLIM
#This file will take the text files and make them into numpy arrays that are
#saved
for curr_dataset in os.listdir(dataDir + enclosing_dir):
    if curr_dataset.split('-')[0] != "mu1e":
        continue
    currDir = dataDir + enclosing_dir + curr_dataset + "/"

    outDir = currDir + "/analysis/"
    #make a directory to save the analysis in
    if not os.path.exists(outDir):
        os.mkdir(outDir)

    #make a dictionary for the timepoint labels
    timepoint_df = pd.read_csv(currDir + 'timepoint_info.tsv', sep = ' ',
                    header= None, names = ['name', 'generation'], index_col = False)

    #a dataframe to save all of the timepoints in
    all_timepoint_genotypes = []

    #Loop through all of the files and get their information.
    for currfile in os.listdir(currDir + "/numpy/"):
        print(currfile, file = sys.stderr)

        #make sure the file has the right type
        if currfile.split('.')[-1] != "npy":
            continue

        #get the timepoint label
        timepoint = currfile.split('_')[-1]
        timepoint = timepoint.split('.')[0]

        #load the array
        coCounts_arr = np.load(currDir + "/numpy/" + currfile)

        #find the segregating sites
        segregatingLoci = zu.find_segregating_diagonal(coCounts_arr, all_seg = True)  
        segregatingLoci['timepoint'] = timepoint_df[timepoint_df['name'] == int(timepoint[-1])]['generation'].tolist()[0]

        #make our dataframe of genotypes
        genotype_df = zu.make_genotype_df(segregatingLoci, coCounts_arr)

        #filter the dataframe 
        genotype_df['timepoint'] = timepoint_df[timepoint_df['name'] == int(timepoint[-1])]['generation'].tolist()[0]
        filtered_genotypes = zu.filter_genotype_df(genotype_df, segregatingLoci, allele_cutoff= CUTOFF, hap_cutoff= SUCCESS)
        all_timepoint_genotypes.append(filtered_genotypes)
        print("Processed one array")


    all_timepoint_genotypes = pd.concat(all_timepoint_genotypes)
    all_timepoint_genotypes.to_pickle(outDir + 'FilteredGenotypes')
    print(all_timepoint_genotypes, file = sys.stderr)



