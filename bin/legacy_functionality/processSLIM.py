import sys
import os
import numpy as np

#For running on cluster
enclosing_dir = '2022_01_25/'
dataDir = '/net/feder/vol1/project/hiv_recombination/data/simulated'

#make a directory for saving output


#Today I am going to be using the data that Alison simulated using SLIM
#This file will take the text files and make them into numpy arrays that are
#saved
for curr_dataset in os.listdir(dataDir + enclosing_dir):
    if curr_dataset.split('-')[0] != "mu1e":
        continue
    currDir = dataDir + enclosing_dir + curr_dataset + "/"
    outDir = currDir + "numpy/"
    if not os.path.exists(outDir):
        os.mkdir(outDir)

    print(outDir)

    #Loop through all of the files and get their information.
    for currfile in os.listdir(currDir):
        if not os.path.isfile(currDir + currfile):
            continue
        if currfile.split('_')[1] != 'formatted':
            continue
        
        #get the timepoint
        time = currfile.split('_')[-1]
        time = time.split('.')[0]

        #load the data into a dataframe
        curr_sim_arr = np.loadtxt(currDir + currfile, dtype = int)
        L = int(np.sqrt(curr_sim_arr.shape[0]/ 36))
        curr_sim_mat = curr_sim_arr.reshape((6, 6, L, L))

        np.save(outDir + currfile.split('.')[0], curr_sim_mat)

    