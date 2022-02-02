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
import zaniniUtil as zu
import os

# # #directories for cluster run
# dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/analysis/neher/'
# dayDir =  '/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/analysis/'
# outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/zanini/researchReports/'

# for running on desktop
dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini/analysis/neher/'
dayDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini/analysis/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/zanini/researchReports/'

#Today I am going to write code to save the reconstitute the results I saved from my Neher analysis and plot them.
#This is similar to the 10-06-2021 file, but I want to add the days between timepoints as the x axis
#I am also adding error bars to my plots based on neher and leitner's approximation

#Next we need to set some things for the run
fragment_list = ['F1','F2', 'F3', 'F4', 'F5', 'F6']
par_list = ['p1', 'p2','p3', 'p4', 'p5', 'p6', 'p7', 'p8', 'p9', 'p10', 'p11']


#viral load bins
NUM_VL_GROUPS = 2
VL_MIN_BIN = 3
VL_BIN_WIDTH = 1.25

#using the custom distance bins

RUNNAME = "viralLoadsSeparate_F1-4"

#We can also load our dataframe for converting days to timepoints
dayToTime = pd.read_csv(dayDir + "DaysToTimepoints.csv")

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


#take the log of viral load
mutation_df['Average_vl'] = np.log10(mutation_df['Average_vl'])
recombination_df['Average_vl'] = np.log10(recombination_df['Average_vl'])

# #label each entry with the distance between timepoints.
# time_diffs_mut = []
# for index, cur_row in mutation_df.iterrows():
#     curr_par = cur_row['Participant']
#     second_time = cur_row['Curr_Timepoint']
#     first_time = cur_row['Last_Timepoint']
#     day_1 = dayToTime[dayToTime['Participant'] == curr_par]
#     day_1 = dayToTime[dayToTime[' Timepoint'] == first_time]
#     day_1 = day_1[' Day'].tolist()[0]
#     day_2 = dayToTime[dayToTime['Participant'] == curr_par]
#     day_2 = dayToTime[dayToTime[' Timepoint'] == second_time]
#     day_2 = day_2[' Day'].tolist()[0]
#     time_diffs_mut.append(day_2 - day_1)

# time_diffs_rec = []
# for index, cur_row in recombination_df.iterrows():
#     curr_par = cur_row['Participant']
#     second_time = cur_row['Curr_Timepoint']
#     first_time = cur_row['Last_Timepoint']
#     day_1 = dayToTime[dayToTime['Participant'] == curr_par]
#     day_1 = dayToTime[dayToTime[' Timepoint'] == first_time]
#     day_1 = day_1[' Day'].tolist()[0]
#     day_2 = dayToTime[dayToTime['Participant'] == curr_par]
#     day_2 = dayToTime[dayToTime[' Timepoint'] == second_time]
#     day_2 = day_2[' Day'].tolist()[0]
#     time_diffs_rec.append(day_2 - day_1)

# mutation_df['Time_Diff'] = time_diffs_mut
# recombination_df['Time_Diff'] = time_diffs_rec
# mutation_df['Dist_x_Time'] = mutation_df['dist'] * mutation_df['Time_Diff']
# recombination_df['Dist_x_Time'] = recombination_df['dist'] * recombination_df['Time_Diff']

# #filter out rows with infinite viral loads (aka not matched to a timepoint)
# mutation_df = mutation_df[mutation_df['Average_vl'] != float('inf')]
# recombination_df = recombination_df[recombination_df['Average_vl'] != float('inf')]

# #make a histogram of the distances multiplied by times
# sns.set(rc={'figure.figsize':(20,5)})
# myplot = sns.FacetGrid(mutation_df, col="Fragment")
# myplot.map_dataframe(sns.histplot, x = 'Dist_x_Time')
# plt.xlim((0, 500000))
# plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
# plt.tight_layout()
# plt.savefig(outDir + "mutation_time_diffHist" + RUNNAME + ".jpg")
# plt.close()

# sns.set(rc={'figure.figsize':(20,5)})
# myplot = sns.FacetGrid(recombination_df, col="Fragment")
# myplot.map_dataframe(sns.histplot, x = 'Dist_x_Time')
# plt.xlim((0, 500000))
# plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
# plt.tight_layout()
# plt.savefig(outDir + "recombination_time_diffHist" + RUNNAME + ".jpg")
# plt.close()





#make a place to store all of the test results
all_frequencies_patients = []

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
            # all_frequencies['Participant'] = curr_par
            # all_frequencies['Fragment'] = curr_fragment
            all_frequencies['Mutation Tests'] = mut_tests
            all_frequencies['Recombination Tests'] = recomb_tests
            all_frequencies['Mutation Successes'] = mut_successes
            all_frequencies['Recombination Successes'] = recomb_successes
            # all_frequencies['Current Timepoint'] = curr_time
            all_frequencies['Avg Viral Load'] = vl_bin_start
            all_frequencies['Fragment'] = frag
            all_frequencies_patients.append(all_frequencies)


#now we can make a dataframe with all our results to plot
all_frequencies_patients = pd.concat(all_frequencies_patients)
all_frequencies_patients['Mut Error'] = 1 / np.sqrt(all_frequencies_patients['Mutation Tests'])
all_frequencies_patients['Recomb Error'] = 1 / np.sqrt(all_frequencies_patients['Recombination Tests'])

for curr_frag in fragment_list:
    frag_data = all_frequencies_patients[all_frequencies_patients['Fragment'] == curr_frag]
    low_vl_data = frag_data[frag_data['Avg Viral Load'] == 3.0]
    high_vl_data = frag_data[frag_data['Avg Viral Load'] == 4.25]

    sns.set(rc={'figure.figsize':(20,5)}, font_scale = 2)
    plt.errorbar(x = low_vl_data['window'], y = low_vl_data['mut_frequencies'], yerr = low_vl_data['Mut Error'], xerr = None, ls = 'none', ecolor = 'gray')
    sns.scatterplot(x = 'window', y = 'mut_frequencies', data = low_vl_data, color = 'gray', label = '3.0 Mutation')
    
    plt.errorbar(x = high_vl_data['window'], y = high_vl_data['mut_frequencies'], yerr = high_vl_data['Mut Error'], xerr = None, ls = 'none', ecolor = 'darkslategray')
    sns.scatterplot(x = 'window', y = 'mut_frequencies', data = high_vl_data, color = 'darkslategray', label = '4.25 Mutation')

    plt.errorbar(x = low_vl_data['window'], y = low_vl_data['recomb_frequencies'], yerr = low_vl_data['Recomb Error'], xerr = None, ls = 'none', ecolor = 'Red')
    sns.scatterplot(x = 'window', y = 'recomb_frequencies', data = low_vl_data, color = 'Red', label = '3.0 Recombination')
    
    plt.errorbar(x = high_vl_data['window'], y = high_vl_data['recomb_frequencies'], yerr = high_vl_data['Recomb Error'], xerr = None, ls = 'none', ecolor = 'Cyan')
    sns.scatterplot(x = 'window', y = 'recomb_frequencies', data = high_vl_data, color = 'Cyan', label = '4.25 Recombination')


    plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
    plt.ylim(-0.1,1.1)
    plt.xlabel("Distance Between Loci X Days Between Timepoints")
    plt.ylabel("Frequency")
    plt.tight_layout()
    plt.savefig(outDir + "mutation_tests" + RUNNAME + curr_frag + ".jpg")
    plt.close()

    # #plot the results for all our participants
    # 
    # sns.set(rc={'figure.figsize':(20,5)})
    # myplot = sns.FacetGrid(all_frequencies_patients, col="Fragment")
    # myplot.map_dataframe(plt.errorbar, x = 'window', y = 'recomb_frequencies', yerr = 'Recomb Error', xerr = None, ls = 'none', ecolor = 'gray')
    # myplot.map_dataframe(sns.scatterplot, x = 'window', y = 'recomb_frequencies', hue = 'Avg Viral Load', palette = 'tab10')    

# plot_groups = np.unique(all_frequencies_patients['Fragment']) 
# #define a color palette index based on column 'B'
# all_frequencies_patients['cind'] = pd.Categorical(all_frequencies_patients['Avg Viral Load'])
# print(all_frequencies_patients['cind'])
# #get the seaborn colour palette and convert to array
# cp = sns.color_palette("rocket", as_cmap = True)
# my_colors = sns.color_palette("rocket", as_cmap = True)
# all_frequencies_patients['cind'] = all_frequencies_patients['Avg Viral Load']/np.linalg.norm(np.unique(all_frequencies_patients['Avg Viral Load']))
# #draw a subplot for each category in column "A"
# fig, axs = plt.subplots(nrows=1, ncols=len(plot_groups), sharey=True)
# for i,ax in enumerate(axs):
#     for vi
#     curr_data = all_frequencies_patients[all_frequencies_patients['Fragment'] == plot_groups[i]]
#     col = list(map (my_colors, curr_data['cind']))
#     print(col)
#     ax.scatter(curr_data['window'], curr_data['recomb_frequencies'], c=col)
#     eb = ax.errorbar(curr_data['window'], curr_data['recomb_frequencies'], yerr=curr_data['Recomb Error'], color = col)

plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
plt.ylim(-0.1,1.1)
plt.xlabel("Distance x Time [BP X Generation]")
plt.ylabel("Frequency")
plt.tight_layout()
plt.savefig(outDir + "recombination_tests" + RUNNAME + ".jpg")
plt.close()


#make a histogram of the average viral loads
sns.set(rc={'figure.figsize':(20,5)})
myplot = sns.FacetGrid(mutation_df, col="Fragment")
myplot.map_dataframe(sns.histplot, x = 'Average_vl')
plt.tight_layout()
plt.savefig(outDir + "mutation_vlHist" + RUNNAME + ".jpg")
plt.close()

#make a histogram of the average viral loads
sns.set(rc={'figure.figsize':(20,5)})
myplot = sns.FacetGrid(recombination_df, col="Fragment")
myplot.map_dataframe(sns.histplot, x = 'Average_vl')
plt.tight_layout()
plt.savefig(outDir + "recombination_vlHist" + RUNNAME + ".jpg")
plt.close()