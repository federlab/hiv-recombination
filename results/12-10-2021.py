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
import os
from scipy import optimize

# #directories for cluster run
# dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/analysis/neher_downsampled_filtered_labeled/'
# outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/zanini/researchReports/bootstrapped_filtered/'

# for running on desktop
dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini/analysis/neher/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/zanini/researchReports/largeBins/'


#Today I just want to try and move functionality to the bin
#This code is to test out whether I have been successful

