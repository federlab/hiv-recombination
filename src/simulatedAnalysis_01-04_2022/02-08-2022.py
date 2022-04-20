import sys
# #for running on cluster
# sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
#for running on desktop
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import os
import numpy as np
import pandas as pd
from statistics import median

#With this script I just want to get the average estimated number of templates
dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini/estimated_depths'

allTemplateNums = []

for currfile in os.listdir(dataDir):
    print(currfile)
    if currfile.split('.')[1] != "tsv":
        continue

    depth_df = pd.read_csv(dataDir +  "/" +currfile, sep = '\t')
    allTemplateNums.extend(depth_df['templates approx'].tolist())
print(allTemplateNums)
print(np.median(allTemplateNums))
print(np.mean(allTemplateNums))