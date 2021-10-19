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

#In this file I am going to conduct the neher analysis, but group the datapoints by viral load.
#To do this we will use the average viral load between the timepoints at which the test
#was conducted.

#directories for cluster run
dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/analysis/'
viralLoadDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/'
outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/zanini/neher_analysis/10-05-2021/Binned_by_VL'
savedData = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/analysis/neher/'

# #for running on desktop
# dataDir = '/Volumes/feder_vol1/home/evromero/2021_hiv-rec/data/zanini/analysis/'
# viralLoadDir = '/Volumes/feder_vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/'
# outDir = '/Volumes/feder_vol1/home/evromero/2021_hiv-rec/results/zanini/neher_analysis/10-05-2021/Binned_by_VL/'
# savedData = '/Volumes/feder_vol1/home/evromero/2021_hiv-rec/data/zanini/analysis/neher/'

#Next we need to set some things for the run
fragment_list = ['F1','F2', 'F3', 'F4', 'F5','F6']
par_list = ['p1', 'p2', 'p3', 'p4', 'p5', 'p6', 'p7', 'p8', 'p9', 'p10', 'p11']
available_files_hap = os.listdir(dataDir + "haplotype_dfs/")
available_files_seg = os.listdir(dataDir + "segregating_Loci/")

# par_list = ['p3', 'p4', 'p5', 'p6', 'p7', 'p8', 'p9', 'p10', 'p11']

NUM_VL_GROUPS = 5
BINWIDTH = 100
MIN_BIN = 0
MAX_BIN = 700
CUTOFF = 0.03
#The frequency a haplotype must have to trigger a test success
SUCCESSFILTER = 0.01
RUNNAME = str(CUTOFF) +  "filtering_success_"  + str(SUCCESSFILTER) +"_binnedVL.png"

#We can start by getting the viral load data
viralLoadData = zu.make_viral_load_df(viralLoadDir)
labeledLoads = []

#This data is in days, but we need to label the timepoints sequentially
for curr_par in par_list:
    curr_data = viralLoadData[viralLoadData['Participant'] == curr_par]
    #pandas gives a warning if you don't explicitely make a copy
    dfWithLabel = curr_data.copy(deep = True)
    dfWithLabel['Timepoint'] = range(1, len(curr_data) + 1)
    labeledLoads.append(dfWithLabel)
viralLoadData = pd.concat(labeledLoads)

#make a dataframe to store our data to plot
all_frequencies_patients = []

#loop through the patients and get their results from the neher analysis
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
        recombination_df, mutation_df = neher.run_analysis(haplotype_df, verbose = False, success_filt = SUCCESSFILTER)

        #label the distances between loci
        recombination_df['dist'] = recombination_df['Locus_2'] - recombination_df['Locus_1']
        mutation_df['dist'] = mutation_df['Locus_2'] - mutation_df['Locus_1']

        #label the test type
        mutation_df['Test_Type'] = 'mutation'
        recombination_df['Test_Type'] = 'recombination'

        #Now get the average viral load between the two timepoints at which the tests were conducted
        #There is probably a clever way to do this in pandas, but I can't think of it
        average_vls = []
        participant_vls = viralLoadData[viralLoadData['Participant'] == curr_par]
        #I am just going to loop through the dataframe and perform each calculation one by one
        for index, row in mutation_df.iterrows():
            time_bef = int(row['Last_Timepoint'])
            success_time = int(row['Curr_Timepoint'])

            #get the viral load at the timepoint that triggered the test
            vl_bef = participant_vls[participant_vls['Timepoint'] == time_bef]
            vl_bef = vl_bef['Viral load [virions/ml]'].tolist()[0]

            #get the viral load at the success timepoint
            vl_success = participant_vls[participant_vls['Timepoint'] == success_time]
            vl_success = vl_success['Viral load [virions/ml]'].tolist()[0]
            average_vls.append(sum([vl_success, vl_bef])/2)
        mutation_df['Average_vl'] = average_vls

        #now repeat the same thing for our recombination df to get the average viral loads
        average_vls = []
        for index, row in recombination_df.iterrows():
            time_bef = int(row['Last_Timepoint'])
            success_time = int(row['Curr_Timepoint'])

            #get the viral load at the timepoint that triggered the test
            vl_bef = participant_vls[participant_vls['Timepoint'] == time_bef]
            vl_bef = vl_bef['Viral load [virions/ml]'].tolist()[0]

            #get the viral load at the success timepoint
            vl_success = participant_vls[participant_vls['Timepoint'] == success_time]
            vl_success = vl_success['Viral load [virions/ml]'].tolist()[0]
            average_vls.append(sum([vl_success, vl_bef])/2)
        recombination_df['Average_vl'] = average_vls

        #write our dataframes to files
        recombination_df.to_csv(savedData + "Recombination_"+ curr_par + "_" + curr_fragment + ".csv")
        mutation_df.to_csv(savedData + "Mutation_"+ curr_par + "_" + curr_fragment + ".csv")

#         #continue if there were no tests run
#         if recombination_df.empty or mutation_df.empty:
#             continue

#         #now we need to group our data by viral load so we can do the distance binning
#         viralLoadRange= recombination_df['Average_vl'].unique().tolist() + mutation_df['Average_vl'].unique().tolist()
#         max_vl = int(max(viralLoadRange))
#         min_vl = int(min(viralLoadRange))

#         #now that we understand the range of the viral loads, we need to make groups
#         #5 groups should correspond to an order of magnitude each
#         vl_groups = range(min_vl, max_vl, (NUM_VL_GROUPS))
#         vl_bin_width = (max_vl-min_vl)/NUM_VL_GROUPS
#         for vl_bin_start in vl_groups:
#             vl_bin_end = vl_bin_start + vl_bin_width
#             #get the data in this group
#             vl_subsetted_recomb = recombination_df[recombination_df['Average_vl'].between(vl_bin_start, vl_bin_end)]
#             vl_subsetted_mut = mutation_df[mutation_df['Average_vl'].between(vl_bin_start, vl_bin_end)]


#             #our next step is to bin by distance
#             #we want to separate out timepoints before we do this binning
#             timepointsList = recombination_df['Curr_Timepoint'].unique()
#             for curr_time in timepointsList:
#                 #get just the tests run at this timepoint
#                 time_subsetted_recomb = vl_subsetted_recomb[vl_subsetted_recomb['Curr_Timepoint'] == curr_time]
#                 time_subsetted_mut = vl_subsetted_mut[vl_subsetted_mut['Curr_Timepoint'] == curr_time]

#                 #create our bins
#                 for bin_start in range(MIN_BIN, MAX_BIN, BINWIDTH):
#                     bin_end = bin_start + BINWIDTH

#                     #make a place to store the frequencies of observing events
#                     mut_frequencies = []
#                     mut_successes = []
#                     recomb_frequencies = []
#                     recomb_successes = []
#                     mut_tests = []
#                     recomb_tests = []

#                     #get all of the datapoints in our bin
#                     curr_recomb = time_subsetted_recomb[time_subsetted_recomb['dist'].between(bin_start, bin_end)]
#                     curr_mut = time_subsetted_mut[time_subsetted_mut['dist'].between(bin_start, bin_end)]

#                     #Calculate the frequencies in each bin
#                     if curr_recomb.shape[0] > 0:
#                         recomb_true = curr_recomb[curr_recomb['Test_Passed'] == True]
#                         recomb_frequencies.append(recomb_true.shape[0]/curr_recomb.shape[0])
#                         recomb_successes.append(recomb_true.shape[0])
#                         recomb_tests.append(curr_recomb.shape[0])

#                     else: 
#                         recomb_frequencies.append(0)
#                         recomb_successes.append(0)
#                         recomb_tests.append(0)

#                     if curr_mut.shape[0] > 0:
#                         mut_true = curr_mut[curr_mut['Test_Passed'] == True]
#                         mut_frequencies.append(mut_true.shape[0]/curr_mut.shape[0])
#                         mut_successes.append(mut_true.shape[0])
#                         mut_tests.append(curr_mut.shape[0])
#                     else: 
#                         mut_frequencies.append(0)
#                         mut_successes.append(0)
#                         mut_tests.append(0)

#                     all_frequencies = pd.DataFrame(list(zip(mut_frequencies,recomb_frequencies)),
#                                             columns = ['mut_frequencies', 'recomb_frequencies'])
#                     all_frequencies['window'] = bin_start
#                     all_frequencies['Participant'] = curr_par
#                     all_frequencies['Fragment'] = curr_fragment
#                     all_frequencies['Mutation Tests'] = mut_tests
#                     all_frequencies['Recombination Tests'] = recomb_tests
#                     all_frequencies['Mutation Successes'] = mut_successes
#                     all_frequencies['Recombination Successes'] = recomb_successes
#                     all_frequencies['Current Timepoint'] = curr_time
#                     all_frequencies['Avg Viral Load'] = vl_bin_start
#                     all_frequencies_patients.append(all_frequencies)

# #now we can make a dataframe with all our results to plot
# all_frequencies_patients = pd.concat(all_frequencies_patients)
# print(all_frequencies_patients)

# #now we can plot our results. at first we'll just color by viral load
# #plot the results for all our participants
# myplot = sns.FacetGrid(all_frequencies_patients, col="Participant")
# sns.set(rc={'figure.figsize':(20,5)})
# myplot.map_dataframe(sns.scatterplot, x = 'window', y = 'mut_frequencies', alpha = 0.5, size = 'Mutation Tests', hue = 'Avg Viral Load')
# # myplot.map_dataframe(sns.lineplot, x = 'window', y = 'mut_frequencies', ci = None)
# plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
# plt.ylim(-0.1,1.1)
# plt.xlim(-10, 500)
# plt.xlabel("Distance Between Loci")
# plt.ylabel("Frequency")
# plt.tight_layout()
# plt.savefig(outDir + "mutation_tests" + RUNNAME)
# plt.close()

# #plot the results for all our participants
# myplot = sns.FacetGrid(all_frequencies_patients, col="Participant")
# sns.set(rc={'figure.figsize':(20,5)})
# myplot.map_dataframe(sns.scatterplot, x = 'window', y = 'recomb_frequencies', alpha = 0.5, size = 'Recombination Tests', hue = 'Avg Viral Load')
# # myplot.map_dataframe(sns.lineplot, x = 'window', y = 'recomb_frequencies', ci = None)
# plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
# plt.ylim(-0.1,1.1)
# plt.xlim(-10, 500)
# plt.xlabel("Distance Between Loci")
# plt.ylabel("Frequency")
# plt.tight_layout()
# plt.savefig(outDir + "recombination_tests" + RUNNAME)
# plt.close()








