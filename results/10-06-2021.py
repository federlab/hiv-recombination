import sys
#for cluster run
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
# #for running on desktop
# sys.path.append('/Volumes/feder_vol1/home/evromero/2021_hiv-rec/bin')
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import neher
import zaniniUtil as zu
import os

# #directories for cluster run
dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/analysis/neher/'
outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/zanini/neher_analysis/10-06-2021/'

#for running on desktop
# dataDir = '/Volumes/feder_vol1/home/evromero/2021_hiv-rec/data/zanini/analysis/neher/'
# outDir = '/Volumes/feder_vol1/home/evromero/2021_hiv-rec/results/zanini/neher_analysis/10-06-2021/'

#Today I am going to write code to save the reconstitute the results I saved from my Neher analysis and plot them.

#Next we need to set some things for the run
fragment_list = ['F1','F2', 'F3', 'F4', 'F5','F6']
par_list = ['p1', 'p2','p3', 'p4', 'p5', 'p6', 'p7', 'p8', 'p9', 'p10', 'p11']

#viral load bins
NUM_VL_GROUPS = 6
VL_MIN_BIN = 0
VL_MAX_BIN = 6

#distance bins
BINWIDTH = 100
MIN_BIN = 0
MAX_BIN = 700
CUTOFF = 0.03
SUCCESS = 0.01
RUNNAME = str(CUTOFF) + '_' + str(SUCCESS) +  "plotbyVL_"

#create lists to store all of the results in
rec_dfs = []
mut_dfs = []

#Loop through all of the files and get their information.
for currfile in os.listdir(dataDir):
    #get the participant and the fragment
    curr_par = currfile.split('_')[1]
    print(curr_par)
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
print(len(mut_dfs))
mutation_df = pd.concat(mut_dfs)
recombination_df = pd.concat(rec_dfs)
print(mutation_df['Participant'].unique())
print('Reconstituted the Data Frames')

#take the log of viral load
mutation_df['Average_vl'] = np.log10(mutation_df['Average_vl'])
recombination_df['Average_vl'] = np.log10(recombination_df['Average_vl'])

#make a place to store all of the test results
all_frequencies_patients = []

#Make bins for the viral load groups
vl_bin_width = int((VL_MAX_BIN-VL_MIN_BIN)/NUM_VL_GROUPS)
vl_groups = range(VL_MIN_BIN, VL_MAX_BIN, vl_bin_width)
print("The viral Load Groups are" + str(vl_groups))

for vl_bin_start in vl_groups:
    vl_bin_end = vl_bin_start + vl_bin_width
    #get the data in this group
    vl_subsetted_recomb = recombination_df[recombination_df['Average_vl'].between(vl_bin_start, vl_bin_end)]
    vl_subsetted_mut = mutation_df[mutation_df['Average_vl'].between(vl_bin_start, vl_bin_end)]


    #our next step is to bin by distance
    # #we want to separate out timepoints before we do this binning
    # timepointsList = recombination_df['Curr_Timepoint'].unique()
    # for curr_time in timepointsList:
    #     #get just the tests run at this timepoint
    #     time_subsetted_recomb = vl_subsetted_recomb[vl_subsetted_recomb['Curr_Timepoint'] == curr_time]
    #     time_subsetted_mut = vl_subsetted_mut[vl_subsetted_mut['Curr_Timepoint'] == curr_time]

    #create our bins
    for bin_start in range(MIN_BIN, MAX_BIN, BINWIDTH):
        bin_end = bin_start + BINWIDTH

        #make a place to store the frequencies of observing events
        mut_frequencies = []
        mut_successes = []
        recomb_frequencies = []
        recomb_successes = []
        mut_tests = []
        recomb_tests = []

        #get all of the datapoints in our bin
        curr_recomb = vl_subsetted_recomb[vl_subsetted_recomb['dist'].between(bin_start, bin_end)]
        curr_mut = vl_subsetted_mut[vl_subsetted_mut['dist'].between(bin_start, bin_end)]


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
        # all_frequencies['Participant'] = curr_par
        # all_frequencies['Fragment'] = curr_fragment
        all_frequencies['Mutation Tests'] = mut_tests
        all_frequencies['Recombination Tests'] = recomb_tests
        all_frequencies['Mutation Successes'] = mut_successes
        all_frequencies['Recombination Successes'] = recomb_successes
        # all_frequencies['Current Timepoint'] = curr_time
        all_frequencies['Avg Viral Load'] = vl_bin_start
        all_frequencies_patients.append(all_frequencies)

#now we can make a dataframe with all our results to plot
all_frequencies_patients = pd.concat(all_frequencies_patients)
print(all_frequencies_patients)

#now we can plot our results. at first we'll just color by viral load
#plot the results for all our participants
# myplot = sns.FacetGrid(all_frequencies_patients, col="Participant")
sns.set(rc={'figure.figsize':(20,5)})
sns.scatterplot(x = 'window', y = 'mut_frequencies',  size = 'Mutation Tests', hue = 'Avg Viral Load', data = all_frequencies_patients)
# myplot.map_dataframe(sns.lineplot, x = 'window', y = 'mut_frequencies', ci = None)
plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
plt.ylim(-0.1,1.1)
plt.xlim(-10, 500)
plt.xlabel("Distance Between Loci")
plt.ylabel("Frequency")
plt.tight_layout()
plt.savefig(outDir + "mutation_tests" + RUNNAME + ".jpg")
plt.close()

#plot the results for all our participants
# myplot = sns.FacetGrid(all_frequencies_patients, col="Participant")
sns.set(rc={'figure.figsize':(20,5)})
sns.scatterplot(x = 'window', y = 'recomb_frequencies', size = 'Recombination Tests', hue = 'Avg Viral Load', data = all_frequencies_patients)
# myplot.map_dataframe(sns.lineplot, x = 'window', y = 'recomb_frequencies', ci = None)
plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
plt.ylim(-0.1,1.1)
plt.xlim(-10, 500)
plt.xlabel("Distance Between Loci")
plt.ylabel("Frequency")
plt.tight_layout()
plt.savefig(outDir + "recombination_tests" + RUNNAME + ".jpg")
plt.close()


#make a histogram of the average viral loads
sns.set(rc={'figure.figsize':(20,5)})
myplot = sns.histplot(data = mutation_df, x = 'Average_vl')
plt.ylabel("Average Viral Load")
plt.tight_layout()
plt.savefig(outDir + "mutation_vlHist" + RUNNAME + ".jpg")
plt.close()

#make a histogram of the average viral loads
sns.set(rc={'figure.figsize':(20,5)})
myplot = sns.histplot(data = recombination_df, x = 'Average_vl')
plt.ylabel("Average Viral Load")
plt.tight_layout()
plt.savefig(outDir + "recombination_vlHist" + RUNNAME + ".jpg")
plt.close()