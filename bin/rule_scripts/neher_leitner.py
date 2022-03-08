import sys
#for running on cluster
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
import os
import neher
import numpy as np
import pandas as pd

dataDir = snakemake.output[0]
print(dataDir, file = sys.stderr)
dataDir = dataDir.split('/')[:-2]
dataDir = "/".join(dataDir)

#The frequency alleles need to be at to be counted as mutation successes
CUTOFF = 0.03

#This script takes in the filtered genotypes and conducts + saves the neher analysis
currFile_gen = dataDir + "/analysis/FilteredGenotypes"
currFile_loc = dataDir + "/analysis/FilteredLoci"

filtered_genotypes = pd.read_pickle(currFile_gen)
filtered_loci = pd.read_pickle(currFile_loc)
#need to make timepoints strings for compatibility with neher analysis file
filtered_genotypes['timepoint'] = list(map(str, filtered_genotypes['timepoint']))
filtered_loci['timepoint'] = list(map(str, filtered_genotypes['timepoint']))

recombination_df = neher.run_analysis(filtered_genotypes)
mutation_df = neher.mutation_analysis(filtered_loci, filtered_loci['frag_len'].unique()[0], CUTOFF)

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

outputDir = currFile_gen.split("/")[:-2]
outputDir = "/".join(outputDir)
print(outputDir, file = sys.stderr)
recombination_df.to_pickle(dataDir + "/neher_res/recombination")
mutation_df.to_pickle(dataDir + "/neher_res/mutation")