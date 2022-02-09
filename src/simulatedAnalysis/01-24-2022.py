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


#Today I am continuing work with the Slim Dataset.
#This script takes in the neher results and plots them
DATASET_LIST = ['r1', 'r2', 'r3']

estimate_df = []

for DATASET_NUM in DATASET_LIST:
    #For running on desktop
    dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/'
    dataDir = dataDir + DATASET_NUM + "/analysis/"
    outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/slimDatasets/'
    outDir = outDir + DATASET_NUM + "/"

    #Start by getting the dataframes with the results of the mutation and recombination tests
    mutation_df = pd.read_pickle(dataDir + "mutation")
    recombination_df = pd.read_pickle(dataDir + "recombination")

    #create our bins
    bin_option1 = [(0,250), (250, 500), (500, 750), (750, 1000), (1000, 2000), (2000, 3000), (3000, 4000), 
    (4000, 5000), (5000, 7500), (7500, 10000), (10000, 15000), (15000, 20000), (20000, 30000), (30000, 40000), (40000, 50000) ,(50000, 60000)]
    #bins that neher and leitner approximately used
    bin_option2 = [(0,5000), (5000, 12500), (12500, 22500), (22500, 30000), (30000, 37500), (37500, 45000), (45000, 52500)]
    bin_lists = [bin_option1, bin_option2]

    for i in range(len(bin_lists)):
        custom_bins = bin_lists[i]
        #Get the test results for the relevant bins
        all_frequencies = []


        for currBin in custom_bins:
            bin_start = currBin[0]
            bin_end = currBin[1]
            #added 1/14/22 plot the bins at the center instead of the start
            bin_center = int((currBin[1] - currBin[0])/2) + bin_start

            #make a place to store the frequencies of observing events
            mut_frequencies = []
            mut_successes = []
            recomb_frequencies = []
            recomb_successes = []
            mut_tests = []
            recomb_tests = []

            #get all of the datapoints in our bin
            curr_recomb = recombination_df[recombination_df['Dist_x_Time'].between(bin_start, bin_end)]
            curr_mut = mutation_df[mutation_df['Dist_x_Time'].between(bin_start, bin_end)]


            #Calculate the frequencies in each bin
            if curr_recomb.shape[0] > 0:
                recomb_true = curr_recomb[curr_recomb['Test_Passed'] == True]
                recomb_frequencies.append(recomb_true.shape[0]/curr_recomb.shape[0])
                recomb_successes.append(recomb_true.shape[0])
                recomb_tests.append(curr_recomb.shape[0])

            else: 
                recomb_frequencies.append(0)
                recomb_successes.append(0)
                recomb_tests.append(0)

            if curr_mut.shape[0] > 0:
                mut_true = curr_mut[curr_mut['Test_Passed'] == True]
                mut_frequencies.append(mut_true.shape[0]/curr_mut.shape[0])
                mut_successes.append(mut_true.shape[0])
                mut_tests.append(curr_mut.shape[0])
            else: 
                mut_frequencies.append(0)
                mut_successes.append(0)
                mut_tests.append(0)

            curr_frequencies = pd.DataFrame(list(zip(mut_frequencies,recomb_frequencies)),
                                columns = ['mut_frequencies', 'recomb_frequencies'])
            curr_frequencies['window'] = bin_end
            curr_frequencies['Mutation Tests'] = mut_tests
            curr_frequencies['Recombination Tests'] = recomb_tests
            curr_frequencies['Mutation Successes'] = mut_successes
            curr_frequencies['Recombination Successes'] = recomb_successes
            all_frequencies.append(curr_frequencies)

        all_frequencies = pd.concat(all_frequencies)


        # ########################### Fitting #######################################
        #get the error bars for each set of tests
        all_frequencies['Mut Error'] = 1 / np.sqrt(all_frequencies['Mutation Tests'])
        all_frequencies['Recomb Error'] =  1/ np.sqrt(all_frequencies['Recombination Tests'])

        #perform our fitting
        coeffs, fit_df = neher.run_neher_fit(c0_fixed = False, lower_bounds = [0,0,0],
                                    upper_bounds = [1,1,1], 
                                    initial = [0.13, 0.1804, 0.000114],
                                    test_results = all_frequencies)

        curr_estimate = neher.estimate_recombination_rate(c0 = coeffs[0], c1 = coeffs[1], c2 = coeffs[2])
        estimate_df.append([curr_estimate, i, DATASET_NUM])
        # ########################### Plotting ###########################################

        #Plot our frequencies with fits
        sns.set(rc={'figure.figsize':(20,5)}, font_scale = 2)
        plt.errorbar(x = all_frequencies['window'], y = all_frequencies['mut_frequencies'],
                yerr = all_frequencies['Mut Error'], xerr = None, ls = 'none', ecolor = 'gray')
        sns.lineplot(x = 'window', y = 'mut_frequencies', data = all_frequencies, color = 'gray', label = 'Mutation Tests')    
        plt.errorbar(x = all_frequencies['window'], y = all_frequencies['recomb_frequencies'],
                yerr = all_frequencies['Recomb Error'], xerr = None, ls = 'none', ecolor = 'red')
        sns.lineplot(x = 'window', y = 'recomb_frequencies', data = all_frequencies, color = 'red', label = 'Recombination Tests')  
        sns.lineplot(x = 'x_vals', y = 'fitted_vals_paper', data = fit_df, color = 'blue', label = 'Neher Fit', linewidth = 2, linestyle = '--')
        sns.lineplot(x = 'x_vals', y = 'fitted_vals', data = fit_df, color = 'black', label = 'Our Fit',linewidth = 2)
        plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
        plt.ylim(-0.1,0.6)
        plt.xlabel("Distance x Time [BP X Generation]")
        plt.ylabel("Frequency")
        plt.tight_layout()
        plt.savefig(outDir + DATASET_NUM + "neherResults_" + str(i) +  ".jpg")
        plt.close()
estimate_df = pd.DataFrame(estimate_df, columns = ['Estimate', 'Binning', 'Dataset'])
print(estimate_df)
