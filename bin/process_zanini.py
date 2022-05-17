import sys
import os
import numpy as np
import shutil

#For running on cluster
textDataDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini/snpPairs/"

#make a directory for saving output
outDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/"


#Today I am going to be using the data that Alison simulated using SLIM
#This file will take the text files and make them into numpy arrays that are
#saved
#Loop through all of the files and get their information.
for currfile in os.listdir(textDataDir):
    if not os.path.isfile(textDataDir + currfile):
        continue
    if currfile.split('.')[-1] != 'npy':
        continue
    
    #get the timepoint, participant, and fragment
    print(currfile)
    curr_info = currfile.split('_')
    timepoint = curr_info[3]
    participant = curr_info[1]
    fragment = curr_info[4]
    fragment = fragment.split('.')[0]

    #move the file and rename it
    move_dir = outDir + participant + "_" + fragment + "/"
    new_name = currfile.split('.')[0] + "_t" + timepoint + "." + currfile.split('.')[1]

    shutil.move(textDataDir + currfile, move_dir + new_name)

