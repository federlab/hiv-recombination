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

#This script is for testing my attempts at fixing the mutation line

# #For running on desktop
dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_01_25/mu1e-04_rho0_Ne1000_M1000_rep1'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_01_25/mu1e-04_rho0_Ne1000_M1000_rep1'

dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/r1'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/slimDatasets/r1'

timepoint_dir = dataDir + "/numpy/"

#make a dictionary for the timepoint labels
timepoint_df = pd.read_csv(dataDir + '/timepoint_info.tsv', sep = ' ',
                header= None, names = ['name', 'generation'], index_col = False)

all_seg = []

#Loop through all of the files and get their information.
for currfile in os.listdir(timepoint_dir):

    #make sure the file has the right type
    if currfile.split('.')[-1] != "npy":
        continue

    #get the timepoint label
    timepoint = currfile.split('_')[-1]
    timepoint = timepoint.split('.')[0]

    #load the array
    coCounts_arr = np.load(timepoint_dir  + currfile)
    fragmentLen = coCounts_arr.shape[-1]

    #find the segregating sites
    segregatingLoci = zu.find_segregating_diagonal(coCounts_arr, all_seg = True)  
    segregatingLoci['timepoint'] = timepoint_df[timepoint_df['name'] == int(timepoint[-1])]['generation'].tolist()[0]
    segregatingLoci['frag_len'] = fragmentLen
    all_seg.append(segregatingLoci)

#get the segregating alleles and conduct mutation tests
all_seg = pd.concat(all_seg, ignore_index= True)
mut_suc_df = neher.mutation_analysis(all_seg, fragmentLen, 0.03)

#the bins for our plot
bin_option2 = [(0,5000), (5000, 12500), (12500, 22500), (22500, 30000), (30000, 37500), (37500, 45000), (45000, 52500)]
all_frequencies = []

for currBin in bin_option2:
    bin_start = currBin[0]
    bin_end = currBin[1]
    #added 1/14/22 plot the bins at the center instead of the start
    bin_center = int((currBin[1] - currBin[0])/2) + bin_start

    #make a place to store the frequencies of observing events
    mut_frequencies = []
    mut_successes = []
    mut_tests = []

    #get all of the datapoints in our bin
    curr_mut_suc = mut_suc_df[mut_suc_df['Dist_x_Time'].between(bin_start, bin_end)]


    #Calculate the frequencies in each bin
    if curr_mut_suc.shape[0] > 0:
        mut_true = curr_mut_suc[curr_mut_suc['Test_Passed'] == True]
        mut_frequencies.append(mut_true.shape[0]/curr_mut_suc.shape[0])
        mut_successes.append(mut_true.shape[0])
        mut_tests.append(curr_mut_suc.shape[0])
    else: 
        mut_frequencies.append(0)
        mut_successes.append(0)
        mut_tests.append(0)
    

    curr_frequencies = pd.DataFrame(list(zip(mut_frequencies)),
                        columns = ['mut_frequencies'])
    curr_frequencies['window'] = bin_end
    curr_frequencies['Mutation Tests'] = mut_tests
    curr_frequencies['Mutation Successes'] = mut_successes

    all_frequencies.append(curr_frequencies)

all_frequencies = pd.concat(all_frequencies, ignore_index= True)

#get the error bars for each set of tests
all_frequencies['Mut Error'] = 1 / np.sqrt(all_frequencies['Mutation Tests'])

#Plot our frequencies with fits
sns.set(rc={'figure.figsize':(20,20)}, font_scale = 2)
fig, axs = plt.subplots(2,2)
axs[0,0].errorbar(x = all_frequencies['window'], y = all_frequencies['mut_frequencies'],
        yerr = all_frequencies['Mut Error'], xerr = None, ls = 'none', ecolor = 'gray')
sns.lineplot(x = 'window', y = 'mut_frequencies', data = all_frequencies, color = 'gray', label = 'Mutation Tests', ax = axs[0,0])    
plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
axs[0,0].set_ylim(0,0.6)
# axs[0,0].set_xlim(0, 50000)
axs[0,0].set_xlabel("Distance x Time [BP X Generation]")
axs[0,0].set_ylabel("Frequency")

#make a histogram of tests per bin
print(all_frequencies['window'])
sns.barplot(x = 'window', y = 'Mutation Tests', data = all_frequencies, ax = axs[0,1])
axs[0,1].set_xlabel("Distance x Time [BP X Generation]")
axs[0,1].set_ylabel("Number of tests")


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
    curr_mut_suc = mut_suc_df[mut_suc_df['dist'].between(bin_start, bin_end)]

    #Calculate the frequencies in each bin
    if curr_mut_suc.shape[0] > 0:
        mut_true = curr_mut_suc[curr_mut_suc['Test_Passed'] == True]
        mut_frequencies.append(mut_true.shape[0]/curr_mut_suc.shape[0])
        mut_successes.append(mut_true.shape[0])
        mut_tests.append(curr_mut_suc.shape[0])
    else: 
        mut_frequencies.append(0)
        mut_successes.append(0)
        mut_tests.append(0)

    curr_frequencies = pd.DataFrame(list(zip(mut_frequencies)),
                        columns = ['mut_frequencies'])
    curr_frequencies['window'] = bin_end
    curr_frequencies['Mutation Tests'] = mut_tests
    curr_frequencies['Mutation Successes'] = mut_successes

    all_frequencies.append(curr_frequencies)

all_frequencies = pd.concat(all_frequencies, ignore_index= True)

#get the error bars for each set of tests
all_frequencies['Mut Error'] = 1 / np.sqrt(all_frequencies['Mutation Tests'])

#Plot our frequencies with fits

axs[1,0].errorbar(x = all_frequencies['window'], y = all_frequencies['mut_frequencies'],
        yerr = all_frequencies['Mut Error'], xerr = None, ls = 'none', ecolor = 'gray')
sns.lineplot(x = 'window', y = 'mut_frequencies', data = all_frequencies, color = 'gray', label = 'Mutation Tests', ax = axs[1,0])    
axs[1,0].set_ylim(0,0.6)
axs[1,0].set_xlim(0,500)
axs[1,0].set_xlabel("Distance [BP]")
axs[1,0].set_ylabel("Frequency")

sns.barplot(x = 'window', y = 'Mutation Tests', data = all_frequencies, ax = axs[1,1])
axs[1,1].set_xlabel("Distance[BP]")
axs[1,1].set_ylabel("Number of tests")

plt.tight_layout()
plt.savefig(outDir + "/mut_res_both.jpg")
plt.close()


