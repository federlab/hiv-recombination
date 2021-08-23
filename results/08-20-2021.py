import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import neher

#In this file we are goint to use the haplotypes and segregating loci in 
#the zanini data to perform the neher and leitner analysis.

dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/analysis/'

fragment_list = ['F1','F2', 'F3', 'F4', 'F5','F6']
par_list = ['p1', 'p2', 'p3', 'p4', 'p5', 'p6', 'p7', 'p8', 'p9', 'p10', 'p11']


for curr_par in par_list:
    for curr_fragment in fragment_list:
        #First we need to load in the dataframes for each patient + fragment.
        par_frag = curr_par +  "_" + curr_fragment
        haplotype_df = pd.read_pickle(dataDir + "haplotype_dfs/haplotypes_" + par_frag + ".pkl")
        segregating_Loci = pd.read_pickle(dataDir + "segregating_Loci/segregating_" + par_frag + ".pkl")
        
        #now we can perform neher analysis on this dataframe
        neher.run_analysis(haplotype_df, segregating_Loci)

    #just try one file for now
        break
    break