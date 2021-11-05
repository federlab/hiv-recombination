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
import matplotlib as mpl
import neher
import zaniniUtil as zu
import os

# #directories for cluster run
dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/analysis/neher/'
dayDir =  '/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/analysis/'
savedData = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/analysis/neher/'

# # for running on desktop
# dataDir = '/Volumes/feder_vol1/home/evromero/2021_hiv-rec/data/zanini/analysis/neher/'
# dayDir = '/Volumes/feder_vol1/home/evromero/2021_hiv-rec/data/zanini/analysis/'
# outDir = '/Volumes/feder_vol1/home/evromero/2021_hiv-rec/results/zanini/neher_analysis/10-26-2021/'

#In this file I am just saving the neher analysis but with the timepoints matched up and labeled

#Next we need to set some things for the run
fragment_list = ['F1','F2', 'F3', 'F4', 'F5','F6']
par_list = ['p1', 'p2','p3', 'p4', 'p5', 'p6', 'p7', 'p8', 'p9', 'p10', 'p11']

#distance bins
BINWIDTH = 5000
MIN_BIN = 0
MAX_BIN = 100000
CUTOFF = 0.03
SUCCESS = 0.01
RUNNAME = str(CUTOFF) + '_' + str(SUCCESS) +  "Truncated100k"

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
mutation_df['Average_vl_log'] = np.log10(mutation_df['Average_vl'])
recombination_df['Average_vl_log'] = np.log10(recombination_df['Average_vl'])

#label each entry with the distance between timepoints.
time_diffs_mut = []
for index, cur_row in mutation_df.iterrows():
    curr_par = cur_row['Participant']
    second_time = cur_row['Curr_Timepoint']
    first_time = cur_row['Last_Timepoint']
    day_1 = dayToTime[dayToTime['Participant'] == curr_par]
    day_1 = dayToTime[dayToTime[' Timepoint'] == first_time]
    day_1 = day_1[' Day'].tolist()[0]
    day_2 = dayToTime[dayToTime['Participant'] == curr_par]
    day_2 = dayToTime[dayToTime[' Timepoint'] == second_time]
    day_2 = day_2[' Day'].tolist()[0]
    time_diffs_mut.append(day_2 - day_1)

time_diffs_rec = []
for index, cur_row in recombination_df.iterrows():
    curr_par = cur_row['Participant']
    second_time = cur_row['Curr_Timepoint']
    first_time = cur_row['Last_Timepoint']
    day_1 = dayToTime[dayToTime['Participant'] == curr_par]
    day_1 = dayToTime[dayToTime[' Timepoint'] == first_time]
    day_1 = day_1[' Day'].tolist()[0]
    day_2 = dayToTime[dayToTime['Participant'] == curr_par]
    day_2 = dayToTime[dayToTime[' Timepoint'] == second_time]
    day_2 = day_2[' Day'].tolist()[0]
    time_diffs_rec.append(day_2 - day_1)

mutation_df['Time_Diff'] = time_diffs_mut
recombination_df['Time_Diff'] = time_diffs_rec
mutation_df['Dist_x_Time'] = mutation_df['dist'] * mutation_df['Time_Diff']
recombination_df['Dist_x_Time'] = recombination_df['dist'] * recombination_df['Time_Diff']

#filter out rows with infinite viral loads (aka not matched to a timepoint)
mutation_df = mutation_df[mutation_df['Average_vl'] != float('inf')]
recombination_df = recombination_df[recombination_df['Average_vl'] != float('inf')]

#save the viral loads matched with timepoints
for curr_par in par_list:
    for curr_fragment in fragment_list:
        recombination_sub_df = recombination_df[recombination_df['Participant'] == curr_par]
        recombination_sub_df = recombination_sub_df[recombination_sub_df['Fragment'] == curr_fragment]
        recombination_sub_df.to_csv(savedData + "Recombination_"+ curr_par + "_" + curr_fragment + ".csv")

        mutation_sub_df = mutation_df[mutation_df['Participant'] == curr_par]
        mutation_sub_df = mutation_sub_df[mutation_sub_df['Fragment'] == curr_fragment]
        mutation_sub_df.to_csv(savedData + "Mutation_"+ curr_par + "_" + curr_fragment + ".csv")