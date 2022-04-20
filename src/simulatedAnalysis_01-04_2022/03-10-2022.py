import sys
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

#This script is for testing my new neher leitner rule which uses the updated mutation line process
CUTOFF = 0.03

# #For running on desktop
dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_02_24/mu1e-05_rho1e-04_Ne10000_M800_rep7'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_02_24/mu1e-05_rho1e-04_Ne10000_M800_rep7'

if not os.path.isdir(outDir):
    os.mkdir(outDir)

# ######################### Start of neher leitner rule ##############################

# #This script takes in the filtered genotypes and conducts + saves the neher analysis
# currFile_gen = dataDir + "/analysis/FilteredGenotypes"
# currFile_loc = dataDir + "/analysis/FilteredLoci"

# filtered_genotypes = pd.read_pickle(currFile_gen)
# filtered_loci = pd.read_pickle(currFile_loc)
# print("Read Pickles")

# #need to make timepoints strings for compatibility with neher analysis file
# filtered_genotypes['timepoint'] = list(map(str, filtered_genotypes['timepoint']))
# filtered_loci['timepoint'] = list(map(str, filtered_loci['timepoint']))

# recombination_df = neher.run_analysis(filtered_genotypes)
# # print("Finished recomb analysis")
# mutation_df = neher.mutation_analysis(filtered_loci, filtered_loci['frag_len'].unique()[0], CUTOFF)

# #convert the distance and time columns into what we need
# recombination_df['Curr_Timepoint'] = list(map(int, recombination_df['Curr_Timepoint']))
# recombination_df['Last_Timepoint'] = list(map(int, recombination_df['Last_Timepoint']))
# mutation_df['Curr_Timepoint'] = list(map(int, mutation_df['Curr_Timepoint']))
# mutation_df['Last_Timepoint'] = list(map(int, mutation_df['Last_Timepoint']))

# recombination_df['dist'] = recombination_df['Locus_2'] - recombination_df['Locus_1']

# recombination_df['Dist_x_Time'] = (recombination_df['Curr_Timepoint'] - recombination_df['Last_Timepoint']) * recombination_df['dist']
# mutation_df['Dist_x_Time'] = (mutation_df['Curr_Timepoint'] - mutation_df['Last_Timepoint']) * mutation_df['dist']

# ######################### End of neher leitner rule ##############################
recombination_df = pd.read_pickle(dataDir + "/neher_res/recombination")
mutation_df = pd.read_pickle(dataDir + "/neher_res/mutation")
print(mutation_df['Test_Passed'].unique())


#the bins for our plot
bin_option2 = [(x, x + 5000) for x in range(0, 55000, 2500)]
all_frequencies = []

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
    all_frequencies.append(curr_frequencies)
    print(curr_frequencies)

all_frequencies = pd.concat(all_frequencies, ignore_index= True)

#get the error bars for each set of tests
all_frequencies['Mut Error'] = 1 / np.sqrt(all_frequencies['Mutation Tests'])
all_frequencies['Recomb Error'] = 1 / np.sqrt(all_frequencies['Recombination Tests'])

#Plot our frequencies with fits
# sns.set(rc={'figure.figsize':(20,30)}, font_scale = 2)
fig, axs = plt.subplots(2,2)
axs[0,0].errorbar(x = all_frequencies['window'], y = all_frequencies['mut_frequencies'],
        yerr = all_frequencies['Mut Error'], xerr = None, ls = 'none', ecolor = 'gray')
sns.lineplot(x = 'window', y = 'mut_frequencies', data = all_frequencies, color = 'gray', label = 'Mutation Tests', ax = axs[0,0])    
axs[0,0].errorbar(x = all_frequencies['window'], y = all_frequencies['recomb_frequencies'],
        yerr = all_frequencies['Recomb Error'], xerr = None, ls = 'none', ecolor = 'red')
sns.lineplot(x = 'window', y = 'recomb_frequencies', data = all_frequencies, color = 'red', label = 'Recombination Tests', ax = axs[0,0])
plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
axs[0,0].set_ylim(0,0.6)
# axs[0,0].set_xlim(0, 50000)
axs[0,0].set_xlabel("Distance x Time [BP X Generation]")
axs[0,0].set_ylabel("Frequency")

#make a histogram of tests per bin
# print(all_frequencies['window'])
# sns.barplot(x = 'window', y = 'Recombination Tests', data = all_frequencies, ax = axs[0,2])
# axs[0,2].set_xlabel("Distance x Time [BP X Generation]")
# axs[0,2].set_ylabel("Number of Recombination Tests")

#make a histogram of tests per bin
sns.barplot(x = 'window', y = 'Mutation Tests', data = all_frequencies, ax = axs[0,1])
axs[0,1].set_xlabel("Distance x Time [BP X Generation]")
axs[0,1].set_ylabel("Number of Mutation Tests")

all_frequencies['Recomb Error'] = 1 / np.sqrt(all_frequencies['Recombination Tests'])

#####################Distance Plot##################

#create our bins
bin_set = [(0, 50), (50, 100), (100, 150), (150, 200), (200, 250), (250, 300), (300, 350), (350, 400), (400, 450), (450, 500) ]


all_frequencies = []

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
    all_frequencies.append(curr_frequencies)

all_frequencies = pd.concat(all_frequencies, ignore_index= True)

#get the error bars for each set of tests
all_frequencies['Mut Error'] = 1 / np.sqrt(all_frequencies['Mutation Tests'])
all_frequencies['Recomb Error'] = 1 / np.sqrt(all_frequencies['Recombination Tests'])

#Plot our frequencies with fits

axs[1,0].errorbar(x = all_frequencies['window'], y = all_frequencies['mut_frequencies'],
        yerr = all_frequencies['Mut Error'], xerr = None, ls = 'none', ecolor = 'gray')
sns.lineplot(x = 'window', y = 'mut_frequencies', data = all_frequencies, color = 'gray', label = 'Mutation Tests', ax = axs[1,0])    
# axs[1,0].errorbar(x = all_frequencies['window'], y = all_frequencies['recomb_frequencies'],
#         yerr = all_frequencies['Recomb Error'], xerr = None, ls = 'none', ecolor = 'red')
# sns.lineplot(x = 'window', y = 'recomb_frequencies', data = all_frequencies, color = 'red', label = 'Recombination Tests', ax = axs[1,0]) 
axs[1,0].set_ylim(0,0.6)
axs[1,0].set_xlim(0,500)
axs[1,0].set_xlabel("Distance [BP]")
axs[1,0].set_ylabel("Frequency")

# sns.barplot(x = 'window', y = 'Recombination Tests', data = all_frequencies, ax = axs[1,2])
# axs[1,2].set_xlabel("Distance[BP]")
# axs[1,2].set_ylabel("Number of Recombination Tests")

#make a histogram of tests per bin
sns.barplot(x = 'window', y = 'Mutation Tests', data = all_frequencies, ax = axs[1,1])
axs[1,1].set_xlabel("Distance x Time [BP X Generation]")
axs[1,1].set_ylabel("Number of Mutation Tests")

plt.tight_layout()
plt.savefig(outDir + "/mut_res_both.jpg")
plt.close()
