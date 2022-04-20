import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
# #for running on desktop
# sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import os
import numpy as np
import neher
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import zaniniUtil as zu

#Today I am trying to determine if something is wrong with the mutation line
#The frequency alleles need to be at to be counted as mutation successes
CUTOFF = 0.03

#For running on cluster
dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_02_24/mu1e-05_rho1e-04_Ne10000_M800_rep'
outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_02_24/mu1e-05_rho1e-04_Ne10000_M800_rep'

# #For running on desktop
# dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_02_24/mu1e-05_rho1e-04_Ne10000_M800_rep'
# outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_02_24/mu1e-05_rho1e-04_Ne10000_M800_rep'

if not os.path.isdir(outDir):
    os.mkdir(outDir)


dataset_dist = []
dataset_dist_time = [] 

all_mut_df = []
all_rec_df = []

for curr_rep in range(1, 11): 

    print(curr_rep)
    recombination_df = pd.read_pickle(dataDir + str(curr_rep) + "/neher_res/recombination")
    mutation_df = pd.read_pickle(dataDir + str(curr_rep) + "/neher_res/mutation")

    currFile_gen = dataDir + str(curr_rep) + "/analysis/FilteredGenotypes"
    currFile_loc = dataDir + str(curr_rep)+ "/analysis/FilteredLoci"

    filtered_genotypes = pd.read_pickle(currFile_gen)
    filtered_loci = pd.read_pickle(currFile_loc)

    #need to make timepoints strings for compatibility with neher analysis file
    filtered_genotypes['timepoint'] = list(map(str, filtered_genotypes['timepoint']))
    filtered_loci['timepoint'] = list(map(str, filtered_loci['timepoint']))
    mutation_df = neher.mutation_analysis(filtered_loci, filtered_loci['frag_len'].unique()[0], CUTOFF)
    mutation_df['Dist_x_Time'] = (mutation_df['Curr_Timepoint'] - mutation_df['Last_Timepoint']) * mutation_df['dist']
    all_mut_df.append(mutation_df)
    all_rec_df.append(recombination_df)

    #the bins for our plot
    bin_option2 = [(x, x + 5000) for x in range(0, 100000, 2500)]
    all_frequencies_dt = []

    for currBin in bin_option2:
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
        all_frequencies_dt.append(curr_frequencies)

    all_frequencies_dt = pd.concat(all_frequencies_dt, ignore_index= True)

    #get the error bars for each set of tests
    all_frequencies_dt['Mut Error'] = 1 / np.sqrt(all_frequencies_dt['Mutation Tests'])
    all_frequencies_dt['Recomb Error'] = 1 / np.sqrt(all_frequencies_dt['Recombination Tests'])
    all_frequencies_dt['rep'] = curr_rep

    #####################Distance Plot##################

    #create our bins
    bin_set = [(0, 50), (50, 100), (100, 150), (150, 200), (200, 250), (250, 300), (300, 350), (350, 400), (400, 450), (450, 500) ]


    all_frequencies_d = []

    for currBin in bin_set:
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
        all_frequencies_d.append(curr_frequencies)

    all_frequencies_d = pd.concat(all_frequencies_d, ignore_index= True)

    #get the error bars for each set of tests
    all_frequencies_d['Mut Error'] = 1 / np.sqrt(all_frequencies_d['Mutation Tests'])
    all_frequencies_d['Recomb Error'] = 1 / np.sqrt(all_frequencies_d['Recombination Tests'])
    all_frequencies_d['rep'] = curr_rep

    dataset_dist.append(all_frequencies_d)
    dataset_dist_time.append(all_frequencies_dt)

all_mut_df = pd.concat(all_mut_df, ignore_index= True)
all_rec_df = pd.concat(all_rec_df, ignore_index= True)
dataset_dist = pd.concat(dataset_dist, ignore_index= True)
dataset_dist_time = pd.concat(dataset_dist_time, ignore_index= True)

sns.scatterplot(x = 'dist', y = 'Dist_x_Time', data = all_mut_df, alpha = 0.15, hue = 'Test_Passed')
plt.xlabel('Distance')
plt.ylabel("Distance x Time [BP X Generation]")
plt.savefig(outDir + "mutation_sanity")
plt.close()


sns.lineplot(x = 'window', y = 'mut_frequencies', data = dataset_dist_time, hue = 'rep')    
# sns.lineplot(x = 'window', y = 'recomb_frequencies', data = dataset_dist_time, hue = 'rep', linestyle = '--' )
plt.ylim(0,0.6)
plt.xlabel("Distance x Time [BP X Generation]")
plt.ylabel("Frequency")
plt.tight_layout()
plt.savefig(outDir + "dist_time")
plt.close()


sns.lineplot(x = 'window', y = 'mut_frequencies', data = dataset_dist, hue = 'rep')    
plt.ylim(0,0.6)
plt.xlabel("Distance")
plt.ylabel("Frequency")
plt.tight_layout()
plt.savefig(outDir + "dist")
plt.close()