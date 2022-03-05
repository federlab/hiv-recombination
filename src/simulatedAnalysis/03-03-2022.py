import sys
# #for running on cluster
# sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
#for running on desktop
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import os
import numpy as np
import neher
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import zaniniUtil as zu

#This script is for testing my attempts at fixing the mutation line

# #For running on desktop
dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/slimDatasets/'

# #For running on cluster
# dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/slimDatasets/'
# outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/slimDatasets/'

#A testing file to calculate genotypes
test_genotype_df = '2022_02_08/mu1e-04_rho0.001_Ne10000_M800_rep1/analysis/FilteredGenotypes'
counts_file_1 = '2022_01_25/mu1e-04_rho0_Ne1000_M1000_rep1/numpy/slim_formatted_t1.npy'
counts_file_2 = '2022_01_25/mu1e-04_rho0_Ne1000_M1000_rep1/numpy/slim_formatted_t2.npy'
counts_file_3 = '2022_01_25/mu1e-04_rho0_Ne1000_M1000_rep1/numpy/slim_formatted_t3.npy'
counts_file_4 = '2022_01_25/mu1e-04_rho0_Ne1000_M1000_rep1/numpy/slim_formatted_t4.npy'

# filtered_genotypes = pd.read_pickle(dataDir + test_genotype_df)
# neher.mutation_analysis(filtered_genotypes, 500)

#setback the filtered genotypes don't appear to be in the form that we would like. 
#specifically, it is already in pairs of loci rather than individual loci. 
#I think we either need to have a different underlying data structure supporting the mutation
#tests, or we need to conglomerate data before.

#We can use the co-counts array and segregating sites 
coCounts_arr1 = np.load(dataDir + counts_file_1)
coCounts_arr2 = np.load(dataDir + counts_file_2)
coCounts_arr3 = np.load(dataDir + counts_file_3)
coCounts_arr4 = np.load(dataDir + counts_file_4)

#Maybe, we actually need the segregating loci
segregatingLoci1 = zu.find_segregating_diagonal(coCounts_arr1, all_seg = True)
segregatingLoci1['timepoint'] = 1
segregatingLoci2 = zu.find_segregating_diagonal(coCounts_arr2, all_seg = True)
segregatingLoci2['timepoint'] = 2
segregatingLoci3 = zu.find_segregating_diagonal(coCounts_arr3, all_seg = True)
segregatingLoci3['timepoint'] = 3
segregatingLoci4 = zu.find_segregating_diagonal(coCounts_arr3, all_seg = True)
segregatingLoci4['timepoint'] = 4
allSegs = pd.concat([segregatingLoci1,segregatingLoci2 , segregatingLoci3, segregatingLoci4], ignore_index = True)
print(allSegs)
neher.mutation_analysis(allSegs, 500, 0.03)