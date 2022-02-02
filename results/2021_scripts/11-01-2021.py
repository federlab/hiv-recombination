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
dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/analysis/neher_downsampled_filtered/'
dayDir =  '/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/analysis/'
savedData = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/analysis/neher_downsampled_filtered_labeled/'
viralLoadDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/'

# # for running on desktop
# dataDir = '/Volumes/feder_vol1/home/evromero/2021_hiv-rec/data/zanini/analysis/neher/'
# dayDir = '/Volumes/feder_vol1/home/evromero/2021_hiv-rec/data/zanini/analysis/'
# outDir = '/Volumes/feder_vol1/home/evromero/2021_hiv-rec/results/zanini/neher_analysis/10-26-2021/'

#In this file I am just saving the neher analysis but with the timepoints matched up and labeled

#Next we need to set some things for the run
fragment_list = ['F1','F2', 'F3', 'F4', 'F5','F6']
par_list = ['p1', 'p2','p3','p5', 'p6', 'p8']
#par_list = ['p9', 'p11']


#distance bins
BINWIDTH = 5000
MIN_BIN = 0
MAX_BIN = 100000
CUTOFF = 0.03
SUCCESS = 0.01
RUNNAME = str(CUTOFF) + '_' + str(SUCCESS) +  "Truncated100k"

#We can start by getting the viral load data
viralLoadData = zu.make_viral_load_df(viralLoadDir)
#We can also load our dataframe for converting days to timepoints
dayToTime = pd.read_csv(dayDir + "DaysToTimepoints.csv")
labeledLoads = []

#We need to add a timepoint column to our viral load data
for curr_par in par_list:
    curr_data = viralLoadData[viralLoadData['Participant'] == curr_par]
    curr_days = dayToTime[dayToTime['Participant'] == curr_par]
    labels = []
    #get the timepoint for the day
    for index, row in curr_data.iterrows():
        date = int(row['Days from infection'])
        currLabel = curr_days[curr_days[' Day'] == date]
        # print("****************", file = sys.stderr)
        # print(curr_par, file = sys.stderr)
        # print(currLabel, file = sys.stderr)
        # print(date, file = sys.stderr)

        labels.append(currLabel[' Timepoint'].tolist()[0])
    labeled_current = curr_data.copy()
    labeled_current['Timepoint'] = labels
    labeledLoads.append(labeled_current)
labeledLoads = pd.concat(labeledLoads)
viralLoadData = labeledLoads
print(viralLoadData, file = sys.stderr)

#create lists to store all of the results in
rec_dfs = []
mut_dfs = []

#Loop through all of the files and get their information.
for currfile in os.listdir(dataDir):
    #get the participant and the fragment
    curr_par = currfile.split('_')[1]
    curr_frag = currfile.split('_')[2]
    curr_frag = curr_frag.split('.')[0]
    curr_sample = currfile.split('_')[3]
    curr_sample = curr_sample.split('.')[0]

    #check if we want to analyze the data from this participant
    if curr_par not in par_list or curr_frag not in fragment_list:
        continue
    
    #read the file into our dataframe
    curr_df = pd.read_csv(dataDir + currfile, index_col= 0)
    curr_df['Participant'] = curr_par
    curr_df['Fragment'] = curr_frag
    curr_df['Sample'] = curr_sample

    #check which type of tests were conducted
    if currfile.split('_')[0] == 'Recombination':
        rec_dfs.append(curr_df)
    else: mut_dfs.append(curr_df)
mutation_df_all = pd.concat(mut_dfs)
recombination_df_all = pd.concat(rec_dfs)

print(mutation_df_all, file = sys.stderr)

#Now get the average viral load between the two timepoints at which the tests were conducted
#There is probably a clever way to do this in pandas, but I can't think of it
for curr_par in par_list:
    for curr_frag in fragment_list:
        mutation_df = mutation_df_all[mutation_df_all['Participant'] == curr_par]
        mutation_df = mutation_df[mutation_df['Fragment'] == curr_frag]
        recombination_df = recombination_df_all[recombination_df_all['Participant'] == curr_par]
        recombination_df = recombination_df[recombination_df['Fragment'] == curr_frag]
        average_vls = []
        participant_vls = viralLoadData[viralLoadData['Participant'] == curr_par]

        #I am just going to loop through the dataframe and perform each calculation one by one
        for index, row in mutation_df.iterrows():
            time_bef = int(row['Last_Timepoint'])
            success_time = int(row['Curr_Timepoint'])
            
            #get the viral load at the timepoint that triggered the test
            vl_bef = participant_vls[participant_vls['Timepoint'] == time_bef]

            if vl_bef.empty:
                average_vls.append(float('inf'))
            else:
                vl_bef = vl_bef['Viral load [virions/ml]'].tolist()[0]
            
            #get the viral load at the success timepoint
            vl_success = participant_vls[participant_vls['Timepoint'] == success_time]
            if vl_success.empty:
                average_vls.append(float('inf'))
            else: 
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
            if vl_bef.empty:
                average_vls.append(float('inf'))

            else: vl_bef = vl_bef['Viral load [virions/ml]'].tolist()[0]

            #get the viral load at the success timepoint
            vl_success = participant_vls[participant_vls['Timepoint'] == success_time]
            if vl_success.empty:
                average_vls.append(float('inf'))
            else:
                vl_success = vl_success['Viral load [virions/ml]'].tolist()[0]
                average_vls.append(sum([vl_success, vl_bef])/2)
        recombination_df['Average_vl'] = average_vls

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
        mutation_df['dist'] = mutation_df['Locus_2'] - mutation_df['Locus_1']
        mutation_df['Dist_x_Time'] = mutation_df['dist'] * mutation_df['Time_Diff']
        recombination_df['dist'] = recombination_df['Locus_2'] - recombination_df['Locus_1']
        recombination_df['Dist_x_Time'] = recombination_df['dist'] * recombination_df['Time_Diff']

        #filter out rows with infinite viral loads (aka not matched to a timepoint)
        mutation_df = mutation_df[mutation_df['Average_vl'] != float('inf')]
        recombination_df = recombination_df[recombination_df['Average_vl'] != float('inf')]
        print(mutation_df, file = sys.stderr)


        # recombination_sub_df = recombination_df[recombination_df['Participant'] == curr_par]
        # recombination_sub_df = recombination_sub_df[recombination_sub_df['Fragment'] == curr_fragment]
        recombination_df.to_csv(savedData + "Recombination_"+ curr_par + "_" + curr_frag + ".csv")
        # print(curr_par, file = sys.stderr)
        # print("Only entered the loop", file = sys.stderr)
        # print(mutation_df, file = sys.stderr)
        # print(mutation_df.columns, file = sys.stderr)
        # print(mutation_df['Participant'].unique(), file = sys.stderr)
        # mutation_sub_df = mutation_df[mutation_df['Participant'] == curr_par]
        # print(mutation_sub_df, file = sys.stderr)
        # mutation_sub_df = mutation_sub_df[mutation_sub_df['Fragment'] == curr_fragment]
        # print(mutation_sub_df, file = sys.stderr)
        mutation_df.to_csv(savedData + "Mutation_"+ curr_par + "_" + curr_frag + ".csv")