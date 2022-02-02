import sys
sys.path.append('/net/gs/vol1/home/evromero/2021_hiv-rec/bin')
import os
import r2Analysis as r2
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

#Today we want to make the plots of R^2 versus distance within each participant
vcfDir = '/net/gs/vol1/home/evromero/2021_hiv-rec/results/vcfs/'
samDir = '/net/gs/vol1/home/evromero/2021_hiv-rec/data/strauli/bam/'
outDir = '/net/gs/vol1/home/evromero/2021_hiv-rec/results/r2_plots/'

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
    
#now make a plot for each group of files
for curr_par in participant_files.keys():
    print(curr_par)
    timepoints = []
    print(curr_par)
    print(participant_files[curr_par])


    #loop through the files and get R^2 values
    for curr_file in participant_files[curr_par]:
        # print("starting loop")
        # print(curr_file)
        # break
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
            timepoints.append(curr_df)
    

    timepoints = pd.concat(timepoints)
    timepoints = timepoints.sort_values(by=['date'])

    ##############THRESHOLD = 50 supporting reads ###############
    timepoints= timepoints[timepoints['support']>=50]

    #make a dataframe of all our results
    sns.set(rc={'figure.figsize':(15,5)})
    myplot = sns.scatterplot(x = 'dist', y = 'r2', size = 'support', hue = 'date', data = timepoints)
    myplot.legend(loc='center left', bbox_to_anchor=(1.25, 0.5), ncol=2)
    plt.ylim(-0.1,1.1)
    plt.xlim(-10,250)
    plt.xlabel("Distance Between Loci")
    plt.ylabel("R^2 Value")
    plt.tight_layout()
    plt.savefig(outDir + curr_par + "filter50")
    plt.close()
    


