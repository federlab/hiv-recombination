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

#This script is to check if we see a pattern in distance data
#This is because all of our curves are showing up flat

# #For running on desktop
dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/slimDatasets/'

# #For running on cluster
# dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/slimDatasets/'
# outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/slimDatasets/'

estimate_df = []

enclosing_dir = '2022_02_24/'

# all_datasets_freqs = []

for curr_dataset in os.listdir(dataDir + enclosing_dir):
    print(curr_dataset, file = sys.stderr)
    if curr_dataset.split('-')[0] != "mu1e":
        continue
    #make a directory for output 
    if not os.path.isdir(outDir + enclosing_dir + curr_dataset):
        os.mkdir(outDir + enclosing_dir + curr_dataset)
    currDir = dataDir + enclosing_dir + curr_dataset + "/neher_res/"

    #Start by getting the dataframes with the results of the mutation and recombination tests
    mutation_df = pd.read_pickle(currDir + "mutation")
    recombination_df = pd.read_pickle(currDir + "recombination")

    #create our bins
    bin_set = [(0, 50), (50, 100), (100, 150), (150, 200), (200, 250), (250, 300), (300, 350), (350, 400), (400, 450), (450, 500) ]


    all_frequencies = []


    for currBin in bin_set:
        print(currBin)
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
        curr_recomb = recombination_df[recombination_df['dist'].between(bin_start, bin_end)]
        curr_mut = mutation_df[mutation_df['dist'].between(bin_start, bin_end)]


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
    all_frequencies['Sim_Rho'] = curr_dataset.split("_")[1]
    all_frequencies['rep'] = curr_dataset.split("_")[-1]
    all_frequencies = all_frequencies.reset_index()

    #get the error bars for each set of tests
    all_frequencies['Mut Error'] = 1 / np.sqrt(all_frequencies['Mutation Tests'])
    all_frequencies['Recomb Error'] =  1/ np.sqrt(all_frequencies['Recombination Tests'])

    #Plot our frequencies with fits
    sns.set(rc={'figure.figsize':(20,5)}, font_scale = 2)
    plt.errorbar(x = all_frequencies['window'], y = all_frequencies['mut_frequencies'],
            yerr = all_frequencies['Mut Error'], xerr = None, ls = 'none', ecolor = 'gray')
    sns.lineplot(x = 'window', y = 'mut_frequencies', data = all_frequencies, color = 'gray', label = 'Mutation Tests')    
    plt.errorbar(x = all_frequencies['window'], y = all_frequencies['recomb_frequencies'],
            yerr = all_frequencies['Recomb Error'], xerr = None, ls = 'none', ecolor = 'red')
    sns.lineplot(x = 'window', y = 'recomb_frequencies', data = all_frequencies, color = 'red', label = 'Recombination Tests')  
    # sns.lineplot(x = 'x_vals', y = 'fitted_vals_paper', data = fit_df, color = 'blue', label = 'Neher Fit', linewidth = 2, linestyle = '--')
    # sns.lineplot(x = 'x_vals', y = 'fitted_vals', data = fit_df, color = 'black', label = 'Our Fit',linewidth = 2)
    plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
    plt.ylim(-0.1,0.6)
    plt.xlabel("Distance x Time [BP X Generation]")
    plt.ylabel("Frequency")
    plt.tight_layout()
    plt.savefig(outDir + enclosing_dir + curr_dataset + "/neherResults_nofit_dist.jpg")
    plt.close()

    break



