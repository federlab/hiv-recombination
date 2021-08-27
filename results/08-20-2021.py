import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import neher
import zaniniUtil as zu

#In this file we are goint to use the haplotypes and segregating loci in 
#the zanini data to perform the neher and leitner analysis.

dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/analysis/'
outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/zanini/neher_analysis'

fragment_list = ['F1','F2', 'F3', 'F4', 'F5','F6']
par_list = ['p1', 'p2', 'p3', 'p4', 'p5', 'p6', 'p7', 'p8', 'p9', 'p10', 'p11']

BINWIDTH = 100
MIN_BIN = 0
MAX_BIN = 700

#make dataframes to store our data to plot
all_frequencies_patients = []

for curr_par in par_list:
    for curr_fragment in fragment_list:
        #First we need to load in the dataframes for each patient + fragment.
        par_frag = curr_par +  "_" + curr_fragment
        haplotype_df = pd.read_pickle(dataDir + "haplotype_dfs/haplotypes_" + par_frag + ".pkl")
        segregating_Loci = pd.read_pickle(dataDir + "segregating_Loci/segregating_" + par_frag + ".pkl")
        
        #now we can perform neher analysis on this dataframe
        #first we will filter out genotypes with alleles at less than 3% frequency
        haplotype_df = zu.filter_genotype_df(haplotype_df, segregating_Loci, cutoff = 0.03)
        recombination_df, mutation_df = neher.run_analysis(haplotype_df, segregating_Loci)

        #label the distances between loci
        recombination_df['dist'] = recombination_df['Locus_2'] - recombination_df['Locus_1']
        mutation_df['dist'] = mutation_df['Locus_2'] - mutation_df['Locus_1']

        #make a place to store the frequencies of observing events
        mut_frequencies = []
        recomb_frequencies = []

        #create our bins
        for bin_start in range(MIN_BIN, MAX_BIN, BINWIDTH):
            bin_end = bin_start + BINWIDTH

            #get all of the datapoints in our bin
            curr_recomb = recombination_df[recombination_df['dist'].between(bin_start, bin_end)]
            curr_mut = mutation_df[mutation_df['dist'].between(bin_start, bin_end)]

            #Calculate the frequencies in each bin
            if curr_recomb.shape[0] > 0:
                recomb_true = curr_recomb[curr_recomb['Test_Passed'] == True]
                recomb_frequencies.append(recomb_true.shape[0]/curr_recomb.shape[0])
            else: mut_frequencies.append(0)

            if curr_mut.shape[0] > 0:
                mut_true = curr_mut[curr_mut['Test_Passed'] == True]
                mut_frequencies.append(mut_true.shape[0]/curr_mut.shape[0])
            else: mut_frequencies.append(0)

        all_frequencies = pd.DataFrame(list(zip(mut_frequencies,recomb_frequencies, list(range(MIN_BIN, MAX_BIN, BINWIDTH)))),
                                columns = ['mut_frequencies', 'recomb_frequencies', 'window'])
        all_frequencies['Participant'] = curr_par
        all_frequencies['Fragment'] = curr_fragment

        all_frequencies_patients.append(all_frequencies)


    #now we can make a dataframe with all our results to plot
    all_frequencies_patients = pd.concat(all_frequencies_patients)

    #plot the results for all our participants
    sns.set(rc={'figure.figsize':(15,5)})
    myplot = sns.scatterplot(x = 'window', y = 'mut_frequencies', hue = 'Participant', data = all_frequencies_patients, alpha = 0.5)
    myplot = sns.scatterplot(x = 'window', y = 'recomb_frequencies', hue = 'Participant', data = all_frequencies_patients, alpha = 0.5)
    myplot.legend(loc='center left', bbox_to_anchor=(1.25, 0.5), ncol=2)
    plt.ylim(-0.1,1.1)
    plt.xlim(-10, MAX_BIN)
    plt.xlabel("Distance Between Loci")
    plt.ylabel("Frequency")
    plt.tight_layout()
    plt.savefig(outDir + "allData_participant")
    plt.close()