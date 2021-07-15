import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
import os
import r2Analysis as r2
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import zaniniUtil as zu

#Now that we can do analysis using the diagonal scanning method, we can plot
dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/snpPairs/'
outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/zanini/'

#parameters for our moving median
#how much to advance the window start by
WINSTEP = 5
#the size of the window
WINSIZE = 20


fragment_list = ['F1','F2', 'F3', 'F4', 'F5','F6']
par_list = ['p1', 'p2', 'p3', 'p4', 'p5', 'p6', 'p7', 'p8', 'p9', 'p10', 'p11']

#keep track of the datafiles for each participant
participant_files = {}

#loop through all of the snp count files and group them by participant and fragment
for file in os.listdir(dataDir):
    #make sure the file has the right type
    if file.split('.')[-1] != "npy":
        continue

    #get the participant and the date
    noExtension = file.split('.')[0]

    #get the participant
    curr_participant = noExtension.split('_')[1]
    curr_fragment = noExtension.split('_')[-1]
    par_frag = curr_participant + "_" + curr_fragment

    #if we have seen files for this participant already add its files to entry
    if par_frag in participant_files.keys():
        participant_files[par_frag] += [file]
    else:
        participant_files[par_frag] = [file]

#Loop through our fragments
for curr_fragment in fragment_list:
    
    #make a list to save all our moving average dataframes in
    all_patients_ave = []
    all_patients_points = []

    #loop through the participants
    for curr_par in par_list:
        timepoints = []
        #for each participant loop through their timepoints
        for curr_file in participant_files[curr_par +  "_" + curr_fragment]:
            #load our current counts
            coCounts_arr = np.load(dataDir + curr_file)
            #find the segregating sites
            segregatingLoci = zu.find_segregating_diagonal(coCounts_arr)

            #continue if there are no segregating sites
            if segregatingLoci.empty:
                continue

            # calculate R^2 values for our locus pairs
            r2List, distList, supportList = r2.calculate_R2_pairCounts(
                coCounts_arr, segregatingLoci)
            curr_df = pd.DataFrame(r2List, columns=['r2'])
            curr_df['dist'] = distList
            curr_df['support'] = supportList
            curr_df['date'] = str(curr_file.split("_")[3])
            curr_df['participant'] = curr_par
            timepoints.append(curr_df)

        #put the data from all the timepoints into a nice dataframe
        timepoints = pd.concat(timepoints)
        timepoints = timepoints.sort_values(by=['date'])
        
        #now we can calculate our moving average
        windowStarts = range(0,int(max(timepoints['dist'])), WINSTEP)
        print(windowStarts, file = sys.stderr)
        patient_aves = []

        ##############THRESHOLD = 50 supporting reads ###############
        timepoints= timepoints[timepoints['support']>=50]

        all_patients_points.append(timepoints)

        #subset the dataframe based on distance
        for i in windowStarts:
            winStart = i
            winEnd = i + WINSIZE
            #get all of the datapoints in our window
            curr_window = timepoints[timepoints['dist'].between(winStart, winEnd)]
            if not curr_window.empty:
                ave_r2 = curr_window['r2'].mean()
                center = winStart + (WINSIZE/2)
                patient_aves.append([center, winStart, winEnd, ave_r2, curr_par])
        
        patient_aves = pd.DataFrame(patient_aves, columns = ['center','window_start', 'window_end', 'average', 'participant'])

        #make a dataframe of all our results
        all_patients_ave.append(patient_aves)

        #plot the results for the current participant
        sns.set(rc={'figure.figsize':(15,5)})
        myplot = sns.scatterplot(x = 'dist', y = 'r2', hue = 'date', data = timepoints, alpha = 0.5)
        myplot.legend(loc='center left', bbox_to_anchor=(1.25, 0.5), ncol=2)
        sns.lineplot(x = 'center', y = 'average', data = patient_aves, linewidth = 3)
        plt.ylim(-0.1,1.1)
        plt.xlim(-10,max(timepoints['dist']))
        plt.xlabel("Distance Between Loci")
        plt.ylabel("R^2 Value")
        plt.tight_layout()
        plt.savefig(outDir + curr_par + "_window_" + str(WINSIZE) + curr_fragment)
        plt.close()

    #convert our list of dataframes across all patients
    all_patients_ave = pd.concat(all_patients_ave)
    all_patients_points = pd.concat(all_patients_points)

    #plot the results for all our participants
    sns.set(rc={'figure.figsize':(15,5)})
    myplot = sns.scatterplot(x = 'dist', y = 'r2', hue = 'participant', data = all_patients_points, alpha = 0.5)
    myplot.legend(loc='center left', bbox_to_anchor=(1.25, 0.5), ncol=2)
    sns.lineplot(x = 'center', y = 'average', data = all_patients_ave, linewidth = 3)
    plt.ylim(-0.1,1.1)
    plt.xlim(-10, max(all_patients_points['dist']))
    plt.xlabel("Distance Between Loci")
    plt.ylabel("R^2 Value")
    plt.tight_layout()
    plt.savefig(outDir + "allPatients_window_" + str(WINSIZE) + curr_fragment)
    plt.close()