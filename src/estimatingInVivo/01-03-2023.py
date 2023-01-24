import numpy as np
import os
import shutil

#Today I am just going to unzip the numpy directories in the zanini data directory

#get the proper names of all the simulated directories
participant_nums = ['p1','p2','p3','p4','p5','p6','p7','p8','p9','p10','p11']
fragments = ['F1', 'F2', 'F3', 'F4', 'F5', 'F6']


dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/'

for curr_dir in os.listdir(dataDir):
    if curr_dir[0] == '.':
        continue
    if 'numpy' in os.listdir(dataDir + curr_dir):
        continue

    #unzip the numpy directory
    os.system("tar -xvzf " + dataDir + curr_dir + "/numpy.tar.gz --directory " + dataDir + curr_dir)
    shutil.move(dataDir + curr_dir + "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/" + 
            curr_dir + "/numpy", dataDir + curr_dir + "/numpy")
    os.system("rm -r " + dataDir + curr_dir + "/Volumes")
   