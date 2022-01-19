import sys
import os
import numpy as np


#Today I am going to be using the data that Alison simulated using SLIM
#This file will take the text files and make them into numpy arrays that are
#saved

# for running on desktop
dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/r3/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/r3/numpy/'

#Loop through all of the files and get their information.
for currfile in os.listdir(dataDir):
    if not os.path.isfile(dataDir + currfile):
        continue
    if currfile.split('_')[1] != 'formatted':
        continue
    
    print(currfile)
    #get the timepoint
    time = currfile.split('_')[-1]
    time = time.split('.')[0]

    #load the data into a dataframe
    curr_sim_arr = np.loadtxt(dataDir + currfile, dtype = int)
    L = int(np.sqrt(curr_sim_arr.shape[0]/ 36))
    curr_sim_mat = curr_sim_arr.reshape((L, L, 6, 6))

    #swap the dimensions since I accidentally gave Alison the wrong order
    swap_sim_mat = np.swapaxes(curr_sim_mat, 0, 2)
    swap_sim_mat = np.swapaxes(swap_sim_mat, 1, 3)

    np.save(outDir + currfile.split('.')[0], swap_sim_mat)

    