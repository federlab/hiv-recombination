import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
import os
import r2Analysis as r2
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

#Now that we can do analysis using the diagonal scanning method, we can plot
dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/snpPairs/'
outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/zanini/'

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
    
