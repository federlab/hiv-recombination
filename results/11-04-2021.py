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

# Today I am going to use the custom bin widths and plot the recombination curves
# they will be binned by viral load and faceted by fragment 

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
par_list = ['p1', 'p2','p3', 'p4', 'p5', 'p6', 'p7', 'p8', 'p9', 'p10', 'p11']
fragment_list = ['F5', 'F6']

#values to truncate our fits at
trunc_vals = list(range(5000, 65000, 5000))
trunc_vals = [60000]

#viral load bins
NUM_VL_GROUPS = 2
VL_MIN_BIN = 3
VL_BIN_WIDTH = 1.25

#filtering values to indicate
CUTOFF = 0.03
SUCCESS = 0.01
RUNNAME = str(CUTOFF) + '_' + str(SUCCESS) +  "vl_binned_custom_improved_f5-6"
#downsample larger distances to get better initial fit.
DOWNSAMPLE_CUTOFF = 1000

#create our distance bins
CUSTOM_BINS = [(0,250), (250, 500), (500, 750), (750, 1000), (1000, 2000), (2000, 3000), (3000, 4000), 
(4000, 5000), (5000, 7500), (7500, 10000), (10000, 15000), (15000, 20000), (20000, 30000), (30000, 40000), (40000, 50000) ,(50000, 60000)]
CUSTOM_BINS = [(0, 1000), (1000, 2000), (2000, 3000), (3000,4000), (4000,5000), (5000, 7500), (7500, 10000), (10000, 15000),
 (15000, 20000), (20000, 30000), (30000, 40000), (40000, 50000) ,(50000, 60000)]

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

#take the log of viral load
mutation_df['Average_vl'] = np.log10(mutation_df['Average_vl'])
recombination_df['Average_vl'] = np.log10(recombination_df['Average_vl'])

#convert the units to basepairs x generation
recombination_df['Dist_x_Time'] = 0.5 * recombination_df['Dist_x_Time']
mutation_df['Dist_x_Time'] = 0.5 * mutation_df['Dist_x_Time']

#Make bins for the viral load groups
vl_groups = [VL_MIN_BIN + VL_BIN_WIDTH * x for x in range (0, NUM_VL_GROUPS)]
print("The viral Load Groups are" + str(vl_groups))

for vl_bin_start in vl_groups:
    vl_bin_end = vl_bin_start + VL_BIN_WIDTH
    #get the data in this group
    vl_subsetted_recomb = recombination_df[recombination_df['Average_vl'].between(vl_bin_start, vl_bin_end)]
    vl_subsetted_mut = mutation_df[mutation_df['Average_vl'].between(vl_bin_start, vl_bin_end)]
    
    #subset by fragment
    for frag in fragment_list:
        vl_frag_recomb = vl_subsetted_recomb[vl_subsetted_recomb['Fragment'] == frag]
        vl_frag_mut = vl_subsetted_mut[vl_subsetted_mut['Fragment'] == frag]

    for currBin in CUSTOM_BINS:
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
        curr_recomb = vl_frag_recomb[vl_frag_recomb['Dist_x_Time'].between(bin_start, bin_end)]
        curr_mut = vl_frag_mut[vl_frag_mut['Dist_x_Time'].between(bin_start, bin_end)]


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
        all_frequencies['Avg Viral Load'] = vl_bin_start
        all_frequencies_patients.append(all_frequencies)

#now we can make a dataframe with all our results to plot
all_frequencies_patients = pd.concat(all_frequencies_patients)
#get the error bars for each set of tests
all_frequencies_patients['Mut Error'] = 1 / np.sqrt(all_frequencies_patients['Mutation Tests'])
all_frequencies_patients['Recomb Error'] =  1/ np.sqrt(all_frequencies_patients['Recombination Tests'])

print(all_frequencies_patients.columns)
lowVL = all_frequencies_patients[all_frequencies_patients['Avg Viral Load'] == 3.0]
highVL = all_frequencies_patients[all_frequencies_patients['Avg Viral Load'] == 4.25]
print(lowVL, file = sys.stderr)
print(highVL, file = sys.stderr)

#plot our results
sns.set(rc={'figure.figsize':(20,10)}, font_scale = 2)
plt.errorbar(x = lowVL['window'], y = lowVL['mut_frequencies'], yerr = lowVL['Mut Error'], xerr = None, ls = 'none', ecolor = 'blue')
sns.lineplot(x = lowVL['window'], y = lowVL['mut_frequencies'], color = 'blue', linestyle = 'dashed', label = "Mutation Low VL")   
plt.errorbar(x = highVL['window'], y = highVL['mut_frequencies'], yerr = highVL['Mut Error'], xerr = None, ls = 'none', ecolor = 'red')
sns.lineplot(x = highVL['window'], y = highVL['mut_frequencies'], color = 'red', linestyle = 'dashed', label = "Mutation High VL")  
plt.errorbar(x = lowVL['window'], y = lowVL['recomb_frequencies'], yerr = lowVL['Recomb Error'], xerr = None, ls = 'none', ecolor = 'blue')
sns.lineplot(x = lowVL['window'], y = lowVL['recomb_frequencies'], color = 'blue', label = "Recombination Low VL")    
plt.errorbar(x = highVL['window'], y = highVL['recomb_frequencies'], yerr = highVL['Recomb Error'], xerr = None, ls = 'none', ecolor = 'red')
sns.lineplot(x = highVL['window'], y = highVL['recomb_frequencies'], color = 'red', label = "Recombination High VL")    

plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
plt.ylim(0,0.6)
plt.xlabel("Distance Between Loci X Days Between Timepoints")
plt.ylabel("Frequency")
plt.tight_layout()
plt.savefig(outDir + "vl_groups" + RUNNAME + ".jpg")
plt.close()
