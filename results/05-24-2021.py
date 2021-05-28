import sys
sys.path.append('/net/gs/vol1/home/evromero/2021_hiv-rec/bin')
import os
import utility
import neher
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import r2Analysis as r2

#Today we want to implement the analysis described by neher and Leitner
vcfDir = '/net/gs/vol1/home/evromero/2021_hiv-rec/results/vcfs/'
samDir = '/net/gs/vol1/home/evromero/2021_hiv-rec/data/strauli/bam/'
outDir = '/net/gs/vol1/home/evromero/2021_hiv-rec/results/r2_plots/'

#To start we need to collect the files for each patient
patient_files = utility.getPatientList(vcfDir)


#Now, for each file, we need to get all of the genotypes.
#get our pairs of loci and the corresponding sam reads
for patient in patient_files.keys():
    #We will store all our genotypes in a pandas dataframe
    all_genotypes = []

    #order and label the timepoints
    filenames = patient_files[patient]
    dates = [int(x.split('y')[1]) for x in filenames]
    dates.sort()
    dateDict = {}
    for i in range(len(dates)): 
        dateDict[dates[i]] = i


    for partic_date in filenames:
    #now that we have our VCF files, we need to start reading them in.
        vcfDF = r2.get_VCF_Loci(vcfDir + partic_date + '.vcf')

        # if there are variants
        if vcfDF is not None :
            #get our pairs of loci and the corresponding sam reads
            genDF = r2.get_sam_reads(vcfDF, samDir + partic_date + '.sorted.bam')
            genDF['date'] = partic_date
            genDF['timepoint'] = dateDict[int(partic_date.split('y')[1])]
        
            #add our dataframe to the list of dataframes we will concatenate
            all_genotypes.append(genDF)
    
    #concatenate all the dataframes across our time points
    all_genotypes = pd.concat(all_genotypes)

    #we need to loop through each timepoint
    timepointList = all_genotypes['timepoint'].unique()
    timepointList.sort()
    #last timepoint left out of loop
    neher.plotPoint(timepointList, all_genotypes)
    break



    
    
