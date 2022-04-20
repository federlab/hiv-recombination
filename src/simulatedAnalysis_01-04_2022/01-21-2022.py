import sys
# #for running on cluster
# sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
#for running on desktop
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import os
import neher
import numpy as np
import pandas as pd


#Today I am continuing work with the Slim Dataset.
#This script takes in the filtered genotypes and conducts + saves the neher analysis
DATASET_NUM = 'r1'

#For running on desktop
dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/'
dataDir = dataDir + DATASET_NUM + "/analysis/"

filtered_genotypes = pd.read_pickle(dataDir + 'FilteredGenotypes')
#need to make timepoints strings for compatibility with neher analysis file
filtered_genotypes['timepoint'] = list(map(str, filtered_genotypes['timepoint']))

recombination_df, mutation_df = neher.run_analysis(filtered_genotypes)

#convert the distance and time columns into what we need
recombination_df['Curr_Timepoint'] = list(map(int, recombination_df['Curr_Timepoint']))
recombination_df['Last_Timepoint'] = list(map(int, recombination_df['Last_Timepoint']))
mutation_df['Curr_Timepoint'] = list(map(int, mutation_df['Curr_Timepoint']))
mutation_df['Last_Timepoint'] = list(map(int, mutation_df['Last_Timepoint']))

recombination_df['dist'] = recombination_df['Locus_2'] - recombination_df['Locus_1']
mutation_df['dist'] = mutation_df['Locus_2'] - mutation_df['Locus_1']

recombination_df['Dist_x_Time'] = (recombination_df['Curr_Timepoint'] - recombination_df['Last_Timepoint']) * recombination_df['dist']
mutation_df['Dist_x_Time'] = (mutation_df['Curr_Timepoint'] - mutation_df['Last_Timepoint']) * mutation_df['dist']

print(recombination_df['Dist_x_Time'])

recombination_df.to_pickle(dataDir + "recombination")
mutation_df.to_pickle(dataDir + "mutation")


