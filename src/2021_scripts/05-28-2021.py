import sys
sys.path.append('/net/gs/vol1/home/evromero/2021_hiv-rec/bin')
import os
import r2Analysis as r2
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

#Today we want to make moving average plots of R^2 over distance in each patient

vcfDir = '/net/gs/vol1/home/evromero/2021_hiv-rec/results/vcfs/'
samDir = '/net/gs/vol1/home/evromero/2021_hiv-rec/data/strauli/bam/'
outDir = '/net/gs/vol1/home/evromero/2021_hiv-rec/results/r2_plots_moving_ave/'

#how much to advance the window start by
WINSTEP = 5
#the size of the window
WINSIZES = [20]

par_list_1 = ['participant1', 'participant4', 'participant7', 'participant9', 'participant10']
par_list_2 = ['participant2', 'participant5', 'participant6']

#keep track of the datafiles for each participant
participant_files = {}

#loop through all of the vcf files and group them by participant
for file in os.listdir(vcfDir):
    #get the participant and the date
    parDate = file.split('.')[0]
    #get the participant
    curr_participant = parDate.split('_')[0]

    #if we have seen files for this participant already add its files to entry
    if curr_participant in participant_files.keys():
        participant_files[curr_participant] += [parDate]
    else:
        participant_files[curr_participant] = [parDate]

#make plots for each specific window size
for WINSIZE in WINSIZES:
    #make a list to save all our moving average dataframes in
    all_patients_ave = []
    all_patients_points = []
        
    #now make a plot for each group of files
    for curr_par in participant_files.keys():
        timepoints = []


        #loop through the files and get R^2 values
        for curr_file in participant_files[curr_par]:
            #now that we have our VCF files, we need to start reading them in.
            vcfDF = r2.get_VCF_Loci(vcfDir + curr_file + '.vcf')

            # if there are variants
            if vcfDF is not None :
                #get our pairs of loci and the corresponding sam reads
                genDF = r2.get_sam_reads(vcfDF, samDir + curr_file + '.sorted.bam')

                #now calculate R^2 values for our locus pairs
                r2List, distList, supportList = r2.r2_all_loci(genDF)
                curr_df = pd.DataFrame(r2List, columns=['r2'])
                curr_df['dist'] = distList
                curr_df['support'] = supportList
                curr_df['date'] = curr_file.split('_')[1]
                curr_df['participant'] = curr_par
                timepoints.append(curr_df)
        

        timepoints = pd.concat(timepoints)
        timepoints = timepoints.sort_values(by=['date'])

        #now we can calculate our moving average
        windowStarts = range(0,250, WINSTEP)
        patient_aves = []

        ##############THRESHOLD = 10 supporting reads ###############
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
        myplot = sns.scatterplot(x = 'dist', y = 'r2', size = 'support', hue = 'date', data = timepoints)
        myplot.legend(loc='center left', bbox_to_anchor=(1.25, 0.5), ncol=2)
        sns.lineplot(x = 'center', y = 'average', data = patient_aves)
        plt.ylim(-0.1,1.1)
        plt.xlim(-10,250)
        plt.xlabel("Distance Between Loci")
        plt.ylabel("R^2 Value")
        plt.tight_layout()
        plt.savefig(outDir + curr_par + "_window_" + str(WINSIZE)  +"_filter_50")
        plt.close()

    #convert our list of dataframes across all patients
    all_patients_ave = pd.concat(all_patients_ave)
    all_patients_points = pd.concat(all_patients_points)
    print(all_patients_ave, file = sys.stderr)

    #plot the results for all our participants
    subset1_ave = all_patients_ave[all_patients_ave['participant'].isin(par_list_1)]
    subset1_poi = all_patients_points[all_patients_points['participant'].isin(par_list_1)]
    print(subset1_ave, file = sys.stderr)
    sns.set(rc={'figure.figsize':(15,5)})
    gr = sns.FacetGrid(subset1_poi,row= 'participant', height=3.5, aspect=2, col_wrap=2)
    gr.map_dataframe(sns.scatterplot, x = 'dist', y = 'r2')
    # gr.map_dataframe(sns.lineplot, x = 'center', y = 'average', data = subset2_ave)
    plt.ylim(-0.1,1.1)
    plt.xlim(-10,250)
    plt.xlabel("Distance Between Loci")
    plt.ylabel("Average R^2 Value")
    plt.tight_layout()
    plt.savefig(outDir + "subset1_" + str(WINSIZE)  +"_filter_50")
    plt.close()

    #plot the results for all our participants
    subset2_ave = all_patients_ave[all_patients_ave['participant'].isin(par_list_2)]
    subset2_poi = all_patients_points[all_patients_points['participant'].isin(par_list_2)]
    print(subset2_ave, file = sys.stderr)
    gr = sns.FacetGrid(subset2_poi, row= 'participant', height=3.5, aspect=2, col_wrap=2)
    gr.map_dataframe(sns.scatterplot, x = 'dist', y = 'r2')
    # gr.map_dataframe(sns.lineplot, x = 'center', y = 'average', data = subset2_ave)
    plt.ylim(-0.1,1.1)
    plt.xlim(-10,250)
    plt.xlabel("Distance Between Loci")
    plt.ylabel("Average R^2 Value")
    plt.tight_layout()
    plt.savefig(outDir + "subset2_" + str(WINSIZE)  +"_filter_50")
    plt.close()