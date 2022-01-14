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
import os



#Today I want to get the autocorrelation method up and running to estimate recombination rate

#I need to start by reading in the data and matching it up with the viral loads.
#Along the way, I'll likely move some functionality to the bin as needed.