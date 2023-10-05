import os
import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import pandas as pd
import autocorrelation as autocorr

#In this file I am just calculating the average sequencing depth in the caskey dataset.
caskey_data = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/caskey/'

all_stat_df_caskey = []
total_loci = 0

all_support = []

#Loop through each of the directory files and calculate the D' ratios
for curr_dir in os.listdir(caskey_data):
    if curr_dir[0] == '.':
        continue
    if not os.path.isdir(caskey_data + curr_dir):
        continue

    inFile = caskey_data + curr_dir + '/linkage/r2_and_D'
    

    #first just try and print the array
    rd_arr = pd.read_pickle(inFile)

    rd_arr.dropna(inplace = True)
    rd_arr['support'] = rd_arr['AB_obs'] + rd_arr['Ab_obs'] + \
         rd_arr['aB_obs'] + rd_arr['ab_obs']

    
    all_support.extend(rd_arr['support'].tolist())

print(all_support)
print('Average support: ', sum(all_support)/len(all_support))
