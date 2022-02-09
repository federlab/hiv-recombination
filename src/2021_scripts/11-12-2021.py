import sys
# #for cluster run
# sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
#for running on desktop
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import neher
import os
from scipy import optimize

# #directories for cluster run
# dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/analysis/neher/'
# dayDir =  '/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/analysis/'
# outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/zanini/researchReports/'

# for running on desktop
dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini/analysis/neher/'
dayDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini/analysis/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/zanini/researchReports/'

# This file plots the neher and leitner analysis binned by participant.

######################### Helper Functions ####################################
def neher_leitner(c0, c1, c2, dDeltaT):
    """The function Neher and Leitner fit to determine recombination rate"""
    return (c0 + (c1 * (1 - np.exp(-c2 * dDeltaT))))

def residual(x0, dDeltaT, data, rec_error):
    #Fix C0 
    c0 = 0.13
    c1 = x0[0]
    c2 = x0[1]

    #this is the equation we are using for our fit
    model = neher_leitner(c0, c1, c2, dDeltaT)
    # print("*****************************", file = sys.stderr)
    # print(model, file = sys.stderr)
    resids = data - model
    # print(data, file = sys.stderr)
    # print(resids, file = sys.stderr)
    weighted_resids = resids * (1 + rec_error)
    return weighted_resids
#################################################################################


#Next we need to set some things for the run
#We are only using fragments 1-4
fragment_list = ['F1','F2', 'F3', 'F4']
par_list = ['p1', 'p2','p3', 'p5', 'p7', 'p8', 'p9', 'p11']

#values to truncate our fits at
trunc_vals = list(range(5000, 65000, 5000))
trunc_vals = [60000]

#distance bins
BINWIDTH = 250
MIN_BIN = 0
MAX_BIN = 60000
CUTOFF = 0.03
SUCCESS = 0.01
RUNNAME = str(CUTOFF) + '_' + str(SUCCESS) +  "fragments1-4"
#downsample larger distances to get better initial fit.
DOWNSAMPLE_CUTOFF = 1000

#create lists to store all of the results in
rec_dfs = []
mut_dfs = []

final_data_df = []


#Loop through all of the files and get their information.
for currfile in os.listdir(dataDir):
    #get the participant and the fragment
    curr_par = currfile.split('_')[1]
    curr_frag = currfile.split('_')[2]
    curr_frag = curr_frag.split('.')[0]

    #check if we want to analyze the data from this participant
    if curr_par not in par_list or curr_frag not in fragment_list:
        continue
    
    #read the file into our dataframe
    curr_df = pd.read_csv(dataDir + currfile, index_col= 0)
    curr_df['Participant'] = curr_par
    curr_df['Fragment'] = curr_frag

    #check which type of tests were conducted
    if currfile.split('_')[0] == 'Recombination':
        rec_dfs.append(curr_df)
    else: mut_dfs.append(curr_df)
mutation_df = pd.concat(mut_dfs)
recombination_df = pd.concat(rec_dfs)

#convert the units to basepairs x generation
recombination_df['Dist_x_Time'] = 0.5 * recombination_df['Dist_x_Time']
mutation_df['Dist_x_Time'] = 0.5 * mutation_df['Dist_x_Time']

for curr_par in par_list:
    #make a place to store all of the test results
    all_frequencies_patients = []

    patient_recomb = recombination_df[recombination_df['Participant'] == curr_par]
    patient_mut = mutation_df[mutation_df['Participant'] == curr_par]


    #create our bins
    custom_bins = [(0,250), (250, 500), (500, 750), (750, 1000), (1000, 2000), (2000, 3000), (3000, 4000), 
    (4000, 5000), (5000, 7500), (7500, 10000), (10000, 15000), (15000, 20000), (20000, 30000), (30000, 40000), (40000, 50000) ,(50000, 60000)]

    for currBin in custom_bins:
        bin_start = currBin[0]
        bin_end = currBin[1]

        #make a place to store the frequencies of observing events
        mut_frequencies = []
        mut_successes = []
        recomb_frequencies = []
        recomb_successes = []
        mut_tests = []
        recomb_tests = []

        #get all of the datapoints in our bin
        curr_recomb = patient_recomb[patient_recomb['Dist_x_Time'].between(bin_start, bin_end)]
        curr_mut = patient_mut[patient_mut['Dist_x_Time'].between(bin_start, bin_end)]


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

        all_frequencies = pd.DataFrame(list(zip(mut_frequencies,recomb_frequencies)),
                                columns = ['mut_frequencies', 'recomb_frequencies'])
        all_frequencies['window'] = bin_start
        all_frequencies['Mutation Tests'] = mut_tests
        all_frequencies['Recombination Tests'] = recomb_tests
        all_frequencies['Mutation Successes'] = mut_successes
        all_frequencies['Recombination Successes'] = recomb_successes
        all_frequencies_patients.append(all_frequencies)

    #now we can make a dataframe with all our results to plot
    all_frequencies_patients = pd.concat(all_frequencies_patients)
    #get the error bars for each set of tests
    all_frequencies_patients['Mut Error'] = 1 / np.sqrt(all_frequencies_patients['Mutation Tests'])
    all_frequencies_patients['Recomb Error'] =  1/ np.sqrt(all_frequencies_patients['Recombination Tests'])

    all_frequencies_patients['Participant'] = curr_par
    final_data_df.append(all_frequencies_patients)



final_data_df = pd.concat(final_data_df)

print(final_data_df['Participant'])
########################### Plotting ###########################################

#Plot our frequencies with fits
sns.set(rc={'figure.figsize':(20,10)}, font_scale = 2)
sns.lineplot(x = 'window', y = 'recomb_frequencies', hue = 'Participant', data = final_data_df, label = 'Recombination Tests')  
sns.scatterplot(x = 'window', y = 'mut_frequencies', data = final_data_df, hue = 'Participant', color = 'black', alpha = 0.5)  
for curr_par in par_list:
    curr_par_df = final_data_df[final_data_df['Participant'] == curr_par]
    plt.errorbar(x = curr_par_df['window'], y = curr_par_df['recomb_frequencies'],
            yerr = curr_par_df['Recomb Error'], xerr = None, ls = 'none', ecolor = 'gray')
plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
plt.xlabel("Distance x Time [BP X Generation]")
plt.ylabel("Frequency")
plt.xlim((0, 50000))
plt.tight_layout()
plt.savefig(outDir + "separated by participant_" + RUNNAME + ".jpg")
plt.close()
