import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import neher
import zaniniUtil as zu
import os

#In this file we are goint to use the haplotypes and segregating loci in 
#the zanini data to perform the neher and leitner analysis.

dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/analysis/'
outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/zanini/neher_analysis/'

fragment_list = ['F1','F2', 'F3', 'F4', 'F5','F6']
par_list = ['p1', 'p2', 'p3', 'p4', 'p5', 'p6', 'p7', 'p8', 'p9', 'p10', 'p11']
available_files_hap = os.listdir(dataDir + "haplotype_dfs/")
available_files_seg = os.listdir(dataDir + "segregating_Loci/")

# par_list = ['p5', 'p6']

BINWIDTH = 100
MIN_BIN = 0
MAX_BIN = 700
CUTOFF = 0.03
#Values to try as the filters for counting a success
SUCCESSFILTERS = [0.01, 0.05, 0.1, 0.2]

for success_filt in SUCCESSFILTERS:
    #name to use for the plots
    RUNNAME = str(CUTOFF) +  "filtering_success_"  + str(success_filt) +"_haps.png"

    #make dataframes to store our data to plot
    all_frequencies_patients = []
    #a dataframe with numbers of supporting reads for tests to use to make a histogram
    support_df = []

    for curr_par in par_list:
        for curr_fragment in fragment_list:
            par_frag = curr_par +  "_" + curr_fragment
            haplotype_file = "haplotypes_" + par_frag + ".pkl"
            loci_file = "segregating_" + par_frag + ".pkl"
            
            #check if there are  files for this
            if haplotype_file not in available_files_hap or loci_file not in available_files_seg:
                continue

            #First we need to load in the dataframes for each patient + fragment.
            haplotype_df = pd.read_pickle(dataDir + "haplotype_dfs/" + haplotype_file)
            segregating_Loci = pd.read_pickle(dataDir + "segregating_Loci/" + loci_file)

            #now we can perform neher analysis on this dataframe
            #first we will filter out genotypes with alleles at less than 3% frequency
            haplotype_df = zu.filter_genotype_df(haplotype_df, segregating_Loci, cutoff = CUTOFF)
            if haplotype_df.empty:
                continue
            recombination_df, mutation_df = neher.run_analysis(haplotype_df, verbose = False, success_filt = success_filt)

            #label the distances between loci
            recombination_df['dist'] = recombination_df['Locus_2'] - recombination_df['Locus_1']
            mutation_df['dist'] = mutation_df['Locus_2'] - mutation_df['Locus_1']

            #get the support for all of the tests into a dataframe so we can plot it for troubleshooting
            mutation_df['Test_Type'] = 'mutation'
            recombination_df['Test_Type'] = 'recombination'
            support_df.append(mutation_df[['dist','Test_Passed','Supporting_Reads', 'Test_Type', 'Success_Freq']])
            support_df.append(recombination_df[['dist','Test_Passed','Supporting_Reads', 'Test_Type', 'Success_Freq']])

            #make a place to store the frequencies of observing events
            mut_frequencies = []
            recomb_frequencies = []
            mut_tests = []
            recomb_tests = []

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
                    recomb_tests.append(curr_recomb.shape[0])

                else: 
                    recomb_frequencies.append(0)
                    recomb_tests.append(0)

                if curr_mut.shape[0] > 0:
                    mut_true = curr_mut[curr_mut['Test_Passed'] == True]
                    mut_frequencies.append(mut_true.shape[0]/curr_mut.shape[0])
                    mut_tests.append(curr_mut.shape[0])
                else: 
                    mut_frequencies.append(0)
                    mut_tests.append(0)

            all_frequencies = pd.DataFrame(list(zip(mut_frequencies,recomb_frequencies, list(range(MIN_BIN, MAX_BIN, BINWIDTH)))),
                                    columns = ['mut_frequencies', 'recomb_frequencies', 'window'])
            all_frequencies['Participant'] = curr_par
            all_frequencies['Fragment'] = curr_fragment
            all_frequencies['Mutation Tests'] = mut_tests
            all_frequencies['Recombination Tests'] = recomb_tests
            all_frequencies_patients.append(all_frequencies)

    #now we can make a dataframe with all our results to plot
    all_frequencies_patients = pd.concat(all_frequencies_patients)
    support_df = pd.concat(support_df)

    #for some further trouble shooting, I'm going to make a plot that indicates the frequency of the
    #haplotype that triggered the success
    sns.set(rc={'figure.figsize':(20,5)})
    myplot = sns.FacetGrid(support_df[support_df['Test_Passed'] == True], col="Test_Type")
    myplot.map_dataframe(sns.scatterplot, y = 'Success_Freq', x = 'dist', hue = 'Test_Passed', size = 'Success_Freq', alpha = 0.5)
    plt.xlim(-10, MAX_BIN)
    plt.xlabel("Distance Between Loci")
    plt.ylabel("Frequency of Success Haplotype")
    plt.tight_layout()
    plt.savefig(outDir + "success_scatter_" + RUNNAME )
    plt.close()


    #I want to do some trouble shooting and make a histogram of the number of reads supporting each test.
    sns.set(rc={'figure.figsize':(20,5)})
    myplot = sns.FacetGrid(support_df, col="Test_Type")
    myplot.map_dataframe(sns.scatterplot, y = 'Supporting_Reads', x = 'dist', hue = 'Test_Passed', alpha = 0.5)
    myplot.set(yscale = 'log')
    plt.xlim(-10, MAX_BIN)
    plt.xlabel("Distance Between Loci")
    plt.ylabel("Supporting Reads")
    plt.tight_layout()
    plt.savefig(outDir + "support_scatter_" + RUNNAME )
    plt.close()

    #plot the results for all our participants
    sns.set(rc={'figure.figsize':(15,5)})
    myplot = sns.FacetGrid(all_frequencies_patients, col="Participant")
    myplot.map_dataframe(sns.scatterplot, x = 'window', y = 'mut_frequencies', alpha = 0.5, size = 'Mutation Tests')
    myplot.map_dataframe(sns.lineplot, x = 'window', y = 'mut_frequencies', ci = None)
    plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
    plt.ylim(-0.1,1.1)
    plt.xlim(-10, 500)
    plt.xlabel("Distance Between Loci")
    plt.ylabel("Frequency")
    plt.tight_layout()
    plt.savefig(outDir + "mutation_tests" + RUNNAME)
    plt.close()

    #plot the results for all our participants
    sns.set(rc={'figure.figsize':(15,5)})
    myplot = sns.FacetGrid(all_frequencies_patients, col="Participant")
    myplot.map_dataframe(sns.scatterplot, x = 'window', y = 'recomb_frequencies', alpha = 0.5, size = 'Recombination Tests')
    myplot.map_dataframe(sns.lineplot, x = 'window', y = 'recomb_frequencies', ci = None)
    plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
    plt.ylim(-0.1,1.1)
    plt.xlim(-10, 500)
    plt.xlabel("Distance Between Loci")
    plt.ylabel("Frequency")
    plt.tight_layout()
    plt.savefig(outDir + "recombination_tests" + RUNNAME)
    plt.close()

    # #make a second plot but size dots by count
    # sns.set(rc={'figure.figsize':(15,5)})
    # myplot = sns.scatterplot(x = 'window', y = 'mut_frequencies', hue = 'Participant', data = all_frequencies_patients, alpha = 0.5,
    #                             palette = 'Greys', size = 'Mutation Tests')
    # myplot = sns.scatterplot(x = 'window', y = 'recomb_frequencies', hue = 'Participant', data = all_frequencies_patients, alpha = 0.5,
    #                             size = 'Recombination Tests')
    # myplot = sns.lineplot(x = 'window', y = 'recomb_frequencies', hue = 'Participant', data = all_frequencies_patients, alpha = 0.5,
    #                     ci = None)
    # myplot.legend(loc='center left', bbox_to_anchor=(1.25, 0.5), ncol=2)
    # plt.ylim(-0.1,1.1)
    # plt.xlim(-10, MAX_BIN)
    # plt.xlabel("Distance Between Loci")
    # plt.ylabel("Frequency")
    # plt.tight_layout()
    # plt.savefig(outDir + "neher_scatter" + RUNNAME )
    # plt.close()

