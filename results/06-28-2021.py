import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
import os
import r2Analysis as r2
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

fragStarts = {}
fragStarts['p1'] = {'F1': 1, 'F2' : 1482, 'F3' : 3124, 'F4' : 4389, 
                  'F5' : 5449, 'F6' : 7266}
fragStarts['p2'] = {'F1': 1, 'F2' : 1477, 'F3' : 3119, 'F4' : 4384, 
                  'F5' : 5441, 'F6' : 7280}
fragStarts['p3'] = {'F1': 1, 'F2' : 1470, 'F3' : 3133, 'F4' : 4398, 
                  'F5' : 5458, 'F6' : 7237}
fragStarts['p4'] = {'F1': 1, 'F2' : 1469, 'F3' : 3111, 'F4' : 4376, 
                  'F5' : 5436, 'F6' : 7275}
fragStarts['p5'] = {'F1': 1, 'F2' : 1469, 'F3' : 3111, 'F4' : 4376, 
                  'F5' : 5424, 'F6' : 7218}
fragStarts['p6'] = {'F1': 1, 'F2' : 1462, 'F3' : 3092, 'F4' : 4357, 
                  'F5' : 5417, 'F6' : 7223}
fragStarts['p7'] = {'F1': 1, 'F2' : 1478, 'F3' : 3126, 'F4' : 4391, 
                  'F5' : 5454, 'F6' : 7287}
fragStarts['p8'] = {'F1': 1, 'F2' : 1468, 'F3' : 3110, 'F4' : 4375, 
                  'F5' : 5432, 'F6' : 7322}
fragStarts['p9'] = {'F1': 1, 'F2' : 1461, 'F3' : 3103, 'F4' : 4368, 
                  'F5' : 5425, 'F6' : 7234}
fragStarts['p10'] = {'F1': 1, 'F2' : 1470, 'F3' : 3112, 'F4' : 4377, 
                  'F5' : 5437, 'F6' : 7291}
fragStarts['p11'] = {'F1': 1, 'F2' : 1466, 'F3' : 3118, 'F4' : 4383, 
                  'F5' : 5443, 'F6' : 7213}




#Today we want to make moving average plots of R^2 over distance in each patient

dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/'
snpPairs = dataDir + "snpPairs/" 
snpCalls = dataDir + "snpCalls/"  
outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/zanini/'

p11_timepoints = [209, 332, 572, 1026, 1396, 1750, 2043]

# snp_loci = r2.get_Zanini_SNP_Loci(testSNPs)

# #we want to run the R^2 analysis on the Zanini data
# coCounts_arr = np.load(testfile)
# #we have a matrix of counts of snp pairs
# r2.calculate_R2_pairCounts(coCounts_arr, snp_loci, 'F6')

par_list_1 = ['p1', 'p2', 'p3', 'p4', 'p5', 'p6', 'p7', 'p8', 'p9', 'p10', 'p11']
# par_list_1 = ['p10']
fragment_list = ['F1','F2', 'F3', 'F4', 'F5','F6']
# fragment_list = ['F2']

#keep track of the datafiles for each participant
participant_files = {}

#loop through all of the snp count files and group them by participant and fragment
for file in os.listdir(snpPairs):
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

#next we need to pair the count files with the files we use to call the variants
call_samples = {}
for par in par_list_1:
    samp_list = []
    for file in os.listdir(snpCalls + par + "/"):
        #make sure the file has the right type
        if file.split('.')[-1] != "tsv":
            continue

        #make an ordered list of sample times
        samp_list.append(int(file.split("_")[0]))
    #sort the time points 
    samp_list.sort()
    call_samples[par] = samp_list

#now we have a dictionary where for sample number i, we look up the relevant
#patient and do samp_list[i-1] to get the corresponding date.

#now we can begin making the plots 
#how much to advance the window start by
WINSTEP = 5
#the size of the window
WINSIZE = 20

#get the snp_loci for each participant
loci_df  = []
for curr_par in par_list_1:
    #loop through the files and get R^2 values
    for curr_t in call_samples[curr_par]:
        #we start by reading in our variants
        snp_loci = r2.get_Zanini_SNP_Loci(
            snpCalls + curr_par + "/" + str(curr_t) + "_days.tsv")

        #filter our loci so they are only segregating (allele frequency > 1%)
        snp_loci = snp_loci[snp_loci["num_seg"] > 1]
        snp_loci["sample_t"] = curr_t
        snp_loci["participant"] = curr_par
        loci_df.append(snp_loci)
    
#put all of our segregating loci into one dataframe
loci_df = pd.concat(loci_df)


#make separate plots for each fragment
for curr_fragment in fragment_list:

    #make a list to save all our moving average dataframes in
    all_patients_ave = []
    all_patients_points = []
    for curr_par in par_list_1:
        timepoints = []

        #get the snps for our participant
        par_loci = loci_df[loci_df['participant'] == curr_par]

        for i in range(len(call_samples[curr_par])):
            sample_t = call_samples[curr_par][i]
            sample_i = i+1

            #now get the relevant files and snps
            par_t_loci = par_loci[par_loci['sample_t'] == sample_t]
            print("Date is " + str(sample_t), file = sys.stderr)
            print(par_t_loci, file = sys.stderr)
            rel_file = snpPairs + "cocounts_" + curr_par + "_sample_" + \
                str(sample_i) + "_" + curr_fragment + ".npy"
            print("Corresponding sample number is " + rel_file, file = sys.stderr)
            if os.path.isfile(rel_file):
                coCounts_arr = np.load(rel_file)

                # if there are variants
                if par_t_loci is not None :
                    # calculate R^2 values for our locus pairs
                    r2List, distList, supportList = r2.calculate_R2_pairCounts(
                        coCounts_arr, par_t_loci, fragStarts[curr_par][curr_fragment], verbose = True)
                    curr_df = pd.DataFrame(r2List, columns=['r2'])
                    curr_df['dist'] = distList
                    curr_df['support'] = supportList
                    curr_df['date'] = sample_t
                    curr_df['participant'] = curr_par
                    timepoints.append(curr_df)
        
        #make sure there are segregating loci
        if len(timepoints) == 0:
            continue
        timepoints = pd.concat(timepoints)
        timepoints = timepoints.sort_values(by=['date'])

        #make sure there are pairs with support
        if timepoints.empty:
            continue
        
        #now we can calculate our moving average
        windowStarts = range(0,int(max(timepoints['dist'])), WINSTEP)
        patient_aves = []

        ##############THRESHOLD = 10 supporting reads ###############
        timepoints= timepoints[timepoints['support']>=50]

        #make sure there are pairs with support
        if timepoints.empty:
            continue

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
        sns.lineplot(x = 'center', y = 'average', data = patient_aves)
        plt.ylim(-0.1,1.1)
        plt.xlim(-10,max(timepoints['dist']))
        plt.xlabel("Distance Between Loci")
        plt.ylabel("R^2 Value")
        plt.tight_layout()
        plt.savefig(outDir + curr_par + "_window_" + str(WINSIZE) + "_" + curr_fragment)
        plt.close()

    #convert our list of dataframes across all patients
    all_patients_ave = pd.concat(all_patients_ave)
    all_patients_points = pd.concat(all_patients_points)

    #plot the results for all our participants
    sns.set(rc={'figure.figsize':(15,5)})
    myplot = sns.scatterplot(x = 'dist', y = 'r2', hue = 'participant', data = all_patients_points, alpha = 0.5)
    myplot.legend(loc='center left', bbox_to_anchor=(1.25, 0.5), ncol=2)
    sns.lineplot(x = 'center', y = 'average', data = all_patients_ave)
    plt.ylim(-0.1,1.1)
    plt.xlim(-10, max(all_patients_points['dist']))
    plt.xlabel("Distance Between Loci")
    plt.ylabel("R^2 Value")
    plt.tight_layout()
    plt.savefig(outDir + "allPatients_window_" + str(WINSIZE) + "_" + curr_fragment)
    plt.close()



