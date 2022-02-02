import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
import os
import r2Analysis as r2
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import zaniniUtil as zu

#Our R^2 analysis worked so now we are going to try with D'

#Now that we can do analysis using the diagonal scanning method, we can plot
dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/snpPairs/'
outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/zanini/d_prime_plots/'
dataOut = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/analysis/'

#parameters for our moving median
#how much to advance the window start by
WINSTEP = 5
#the size of the window
WINSIZE = 20
WINSIZE_ALL = 10


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
        par_frag = curr_par +  "_" + curr_fragment
        par_frag_results = []

        #if there isn't a file for this combination of participant and fragment
        if par_frag not in participant_files.keys():
            continue

        #for each participant loop through their timepoints
        for curr_file in participant_files[par_frag]:
            #load our current counts
            coCounts_arr = np.load(dataDir + curr_file)
            #find the segregating sites
            segregatingLoci = zu.find_segregating_diagonal(coCounts_arr)

            #continue if there are no segregating sites
            if segregatingLoci.empty:
                continue

            # calculate D' values for our locus pairs
            DList, distList, supportList, resultsDf = r2.calculate_R2_pairCounts(
                coCounts_arr, segregatingLoci, statistic = 'D', saveData = True )
            #add our results to be saved
            resultsDf['date'] = str(curr_file.split("_")[3])
            par_frag_results.append(resultsDf)
            #calculate the R^2 value just to save
            curr_df = pd.DataFrame(DList, columns=['D'])
            curr_df['dist'] = distList
            curr_df['support'] = supportList
            curr_df['date'] = str(curr_file.split("_")[3])
            curr_df['participant'] = curr_par
            timepoints.append(curr_df)
        #Save our calculations in a csv file
        par_frag_results = pd.concat(par_frag_results)
        par_frag_results.to_csv(dataOut + 'RandD_' + par_frag + ".csv", index = False)


        #put the data from all the timepoints into a nice dataframe
        timepoints = pd.concat(timepoints)
        timepoints = timepoints.sort_values(by=['date'])
        
        #now we can calculate our moving average
        windowStarts = range(0,int(max(timepoints['dist'])), WINSTEP)
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
                ave_D = curr_window['D'].mean()
                center = winStart + (WINSIZE/2)
                patient_aves.append([center, winStart, winEnd, ave_D, curr_par])
        
        patient_aves = pd.DataFrame(patient_aves, columns = ['center','window_start', 'window_end', 'average', 'participant'])

        #make a dataframe of all our results
        all_patients_ave.append(patient_aves)

        #plot the results for the current participant
        sns.set(rc={'figure.figsize':(15,5)})
        myplot = sns.scatterplot(x = 'dist', y = 'D', hue = 'date', data = timepoints, alpha = 0.5)
        myplot.legend(loc='center left', bbox_to_anchor=(1.25, 0.5), ncol=2)
        sns.lineplot(x = 'center', y = 'average', data = patient_aves, linewidth = 3)
        plt.ylim(-0.1,1.1)
        plt.xlim(-10,max(timepoints['dist']))
        plt.xlabel("Distance Between Loci")
        plt.ylabel("D' Value")
        plt.tight_layout()
        plt.savefig(outDir + curr_par + "_window_" + str(WINSIZE) + curr_fragment)
        plt.close()

    #convert our list of dataframes across all patients
    all_patients_ave = pd.concat(all_patients_ave)
    all_patients_points = pd.concat(all_patients_points)

    #plot the results for all our participants
    sns.set(rc={'figure.figsize':(15,5)})
    myplot = sns.scatterplot(x = 'dist', y = 'D', hue = 'participant', data = all_patients_points, alpha = 0.5)
    myplot.legend(loc='center left', bbox_to_anchor=(1.25, 0.5), ncol=2)
    sns.lineplot(x = 'center', y = 'average', data = all_patients_ave, linewidth = 3)
    plt.ylim(-0.1,1.1)
    plt.xlim(-10, max(all_patients_points['dist']))
    plt.xlabel("Distance Between Loci")
    plt.ylabel("D' Value")
    plt.tight_layout()
    plt.savefig(outDir + "allPatients_window_" + str(WINSIZE_ALL) + curr_fragment)
    plt.close()