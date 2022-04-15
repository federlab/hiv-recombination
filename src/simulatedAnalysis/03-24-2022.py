import sys
# sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
#for running on desktop
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import optimize

#Today I am going to take one of the linkage files and try to plot the 
#autocorrelation of the D statistic
#I need to plot the D values, in case we want to do thresholding.
#I think the easiest thing would be a histogram 

#For running on cluster
dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_02_24/'
outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_02_24/'

dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_02_24/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_02_24/'

estimate_df = [] 

for curr_data in os.listdir(dataDir):
    print(curr_data)
    #only get the data directories, not hidden files
    if curr_data[0] == '.':
        continue
    run_info = curr_data.split('_')
    sim_rho = run_info[1]
    sim_rho = sim_rho[3:]
    rep = run_info[-1]

    #make a place to store our output
    currOut = outDir + curr_data
    if not os.path.exists(currOut):
        os.mkdir(currOut)
    
    linkage_file = dataDir + curr_data + "/linkage/r2_and_D"

    #first just try and print the array
    rd_arr = pd.read_pickle(linkage_file)
    rd_arr.dropna(inplace = True)
    rd_arr['support'] = rd_arr['AB_obs'] + rd_arr['Ab_obs'] + rd_arr['aB_obs'] + rd_arr['ab_obs']
    rd_arr['max_freq'] = rd_arr[['p_A','p_B']].max(axis=1)



    rd_arr = rd_arr.loc[(rd_arr['AB_obs'] > 0) & (rd_arr['aB_obs'] > 0)  & (rd_arr['Ab_obs'] > 0)  & (rd_arr['ab_obs'] > 0) ]
    


    sns.histplot(x = 'd_prime', data = rd_arr)
    plt.savefig(currOut + "/d_suss_hist.jpg")
    plt.close()

    sns.scatterplot(x = 'd_prime', y = 'max_freq', data = rd_arr, alpha = 0.05)
    plt.savefig(currOut + "/d_freq_scatter.jpg")
    plt.close()
    break
