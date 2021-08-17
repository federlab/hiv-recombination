import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import zaniniUtil as zu
import neher

#This script is to try and get Neher and Leitner's analysis 
#from 2009 up and running for the data from Zanini

dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/snpPairs/'


#my idea, we can probably do the analysis for each pair of loci at a time. (loop through loci pairs
# in the outer loop rather than looping through time points)
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
    #make a list to save our genotypes in (this will become a dataframe)
    fragment_genotypes = []
    segregatingLoci_all = []

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
            segregatingLoci = zu.find_segregating_diagonal(coCounts_arr, True)

            #filter out errors
            #I think to do this we just need to check that allele 2 has frequency greater than 3%
            segregatingLoci = segregatingLoci[segregatingLoci['freq_2'] > 0.03]
            segregatingLoci.reset_index(drop=True, inplace=True)


            #continue if there are no segregating sites
            if segregatingLoci.empty:
                continue
            genotype_df = zu.make_genotype_df(segregatingLoci, coCounts_arr)
            #we need to label the timepoint here
            timepoint_name = curr_file.split('_')[3]
            genotype_df['timepoint'] = timepoint_name
            segregatingLoci['timepoint'] = timepoint_name
            fragment_genotypes.append(genotype_df)
            segregatingLoci_all.append(segregatingLoci)

        #Now we have a dataframe of all the haplotypes
        fragment_genotypes = pd.concat(fragment_genotypes)
        segregatingLoci_all = pd.concat(segregatingLoci_all)

        #Next we need to perform neher leitner analysis using these dataframes
        neher.plotPoint(fragment_genotypes, segregatingLoci_all )
        break
    break