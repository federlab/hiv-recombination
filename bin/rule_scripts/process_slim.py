import sys
import os
import numpy as np
import shutil

#For running on cluster
textDataDir = snakemake.input[0]
textDataDir = textDataDir.split('/')[:-1]
textDataDir = "/".join(textDataDir)
textDataDir = textDataDir + "/"

#make a directory for saving output
outDir = snakemake.output[0]
outDir = outDir.split('/')[:-1]
outDir = "/".join(outDir)
if not os.path.exists(outDir):
    os.mkdir(outDir)
if not os.path.exists(outDir + "/numpy/"):
    os.mkdir(outDir + "/numpy/")

#copy the timepoint info into that directory
if not os.path.exists(outDir + "/timepoint_info.tsv"):
    shutil.copy(textDataDir + "/timepoint_info.tsv", outDir + "/timepoint_info.tsv")

#Today I am going to be using the data that Alison simulated using SLIM
#This file will take the text files and make them into numpy arrays that are
#saved
#Loop through all of the files and get their information.
for currfile in os.listdir(textDataDir):
    if not os.path.isfile(textDataDir + currfile):
        continue
    if currfile.split('_')[1] != 'formatted':
        continue
    
    #get the timepoint
    time = currfile.split('_')[-1]
    time = time.split('.')[0]

    #load the data into a dataframe
    curr_sim_arr = np.loadtxt(textDataDir + currfile, dtype = int)
    L = int(np.sqrt(curr_sim_arr.shape[0]/ 36))
    curr_sim_mat = curr_sim_arr.reshape((6, 6, L, L))

    np.save(outDir + "/numpy/" + currfile.split('.')[0], curr_sim_mat)

    