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
# dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/analysis/neher_downsampled_filtered_labeled/'
# outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/zanini/researchReports/bootstrapped_filtered/'

# for running on desktop
dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini/analysis/neher/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/zanini/testingBinning/'


#Today I just want to try and move functionality to the bin
#This code is to test out whether I have been successful

#Next we need to set some things for the run
#We are only using fragments 1-4
fragment_list = ['F1', 'F2', 'F3', 'F4']
# fragment_list = ['F5']
par_list = ['p1', 'p2','p3', 'p4', 'p5', 'p6', 'p7', 'p8', 'p9', 'p10', 'p11']

#These parameters aren't used here but were used when processing the data
CUTOFF = 0.03
SUCCESS = 0.01
RUNNAME = str(CUTOFF) + '_' + str(SUCCESS) +  "moving_functionality_test"

#create lists to store all of the results in
rec_dfs = []
mut_dfs = []

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

#make a place to store all of the test results
all_frequencies_patients = []

#convert the units to basepairs x generation
recombination_df['Dist_x_Time'] = 0.5 * recombination_df['Dist_x_Time']
mutation_df['Dist_x_Time'] = 0.5 * mutation_df['Dist_x_Time']

#create our bins
custom_bins = [(0,250), (250, 500), (500, 750), (750, 1000), (1000, 2000), (2000, 3000), (3000, 4000), 
(4000, 5000), (5000, 7500), (7500, 10000), (10000, 15000), (15000, 20000), (20000, 30000), (30000, 40000), (40000, 50000) ,(50000, 60000)]
#bins that neher and leitner approximately used
custom_bins = [(0,5000), (5000, 12500), (12500, 22500), (22500, 30000), (30000, 37500), (37500, 45000), (45000, 52500)]

for currBin in custom_bins:
    bin_start = currBin[0]
    bin_end = currBin[1]
    #added 1/14/22 plot the bins at the center instead of the start
    bin_center = int((currBin[1] - currBin[0])/2) + bin_start
    print(bin_center)

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

    all_frequencies = pd.DataFrame(list(zip(mut_frequencies,recomb_frequencies)),
                            columns = ['mut_frequencies', 'recomb_frequencies'])
    all_frequencies['window'] = bin_end
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

#perform our fitting
coeffs, fit_df = neher.run_neher_fit(c0_fixed = False, lower_bounds = [0,0,0],
                            upper_bounds = [1,1,1], 
                            initial = [0.13, 0.1804, 0.000114],
                            test_results = all_frequencies_patients)

print(coeffs, file = sys.stderr)

########################### Plotting ###########################################

#Plot our frequencies with fits
sns.set(rc={'figure.figsize':(20,5)}, font_scale = 2)
plt.errorbar(x = all_frequencies_patients['window'], y = all_frequencies_patients['mut_frequencies'],
        yerr = all_frequencies_patients['Mut Error'], xerr = None, ls = 'none', ecolor = 'gray')
sns.lineplot(x = 'window', y = 'mut_frequencies', data = all_frequencies_patients, color = 'gray', label = 'Mutation Tests')    
plt.errorbar(x = all_frequencies_patients['window'], y = all_frequencies_patients['recomb_frequencies'],
        yerr = all_frequencies_patients['Recomb Error'], xerr = None, ls = 'none', ecolor = 'red')
sns.lineplot(x = 'window', y = 'recomb_frequencies', data = all_frequencies_patients, color = 'red', label = 'Recombination Tests')  
sns.lineplot(x = 'x_vals', y = 'fitted_vals_paper', data = fit_df, color = 'blue', label = 'Neher Fit', linewidth = 2, linestyle = '--')
sns.lineplot(x = 'x_vals', y = 'fitted_vals', data = fit_df, color = 'black', label = 'Our Fit',linewidth = 2)
plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
plt.ylim(-0.1,0.6)
plt.xlabel("Distance x Time [BP X Generation]")
plt.ylabel("Frequency")
plt.tight_layout()
plt.savefig(outDir + "allTogether_" + RUNNAME + ".jpg")
plt.close()