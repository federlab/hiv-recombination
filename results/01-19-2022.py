import sys
# #for cluster run
# sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
#for running on desktop
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import neher
import zaniniUtil as zu


#for running on desktop
dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/slimDatasets/'

DATASET_NUM = 'r1'

dataDir = dataDir + DATASET_NUM + "/numpy/"
outDir = outDir + DATASET_NUM + "/"

#read in the data
for currfile in os.listdir(dataDir):
    #make sure the file has the right type
    if currfile.split('.')[-1] != "npy":
        continue

    #get the timepoint label
    timepoint = currfile.split('_')[-1]
    timepoint = timepoint.split('.')[0]

    #load the array
    coCounts_arr = np.load(dataDir + currfile)

    #find the segregating sites
    segregatingLoci = zu.find_segregating_diagonal(coCounts_arr)    

