import os
import sys
import numpy as np
import pandas as pd
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import autocorrelation as autocorr

#In this file I am loading in the caskey viral load data and plotting some summary statistics

#For running on desktop
VL_dataDir = '/Volumes/feder-vol1/home/evromero/2023_hiv-bnabs/data/caskey2017/caskey_VL/'
caskey_data = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/caskey/'

#A dictionary to map timepoints to days
timepoint_dict = {'preinf' : 0,
                  'D0' : 0,
                  'W1' : 7 * 0.5,
                  'W2' : 14 * 0.5,
                  'W4' : 28 * 0.5,
                  'W8' : 56 * 0.5,
                  'W12' : 84 * 0.5,
                  'W16' : 112 * 0.5,
                  'W20' : 140 * 0.5,
                  'W24' : 168 * 0.5}

# Make a dictionary of the viral load data
vl_df = pd.read_csv(VL_dataDir + 'caskey_VL_data.tsv', sep='\t')
vl_dict = {}
for index, row in vl_df.iterrows():
    time = row['Timepoint']
    time = timepoint_dict[time]
    indiv = row['Individual']
    time_indiv = str(time) + '_' + indiv
    vl_measure = row['Viral Load (cp/mL)'].replace(',', '')
    vl_dict[time_indiv] = int(vl_measure)

all_stat_df_caskey = []
total_loci = 0

#Loop through each of the directory files and calculate the D' ratios
for curr_dir in os.listdir(caskey_data):
    if curr_dir[0] == '.':
        continue
    if not os.path.isdir(caskey_data + curr_dir):
        continue
    
    inFile = caskey_data + curr_dir + '/linkage/r2_and_D'

    #Get the D' ratios
    stat_df = autocorr.calculate_d_ratios(inFile)

    stat_df = stat_df[stat_df['Dist_X_Time'] <= 50000]
    if len(stat_df) == 0:
        continue
    
    stat_df['Participant'] = curr_dir
    all_stat_df_caskey.append(stat_df)

all_stat_df_caskey = pd.concat(all_stat_df_caskey, ignore_index=True)
all_stat_df_caskey['Time_Par1'] = all_stat_df_caskey['Time_1'].astype(str) + '_' + all_stat_df_caskey['Participant']
all_stat_df_caskey['Time_Par2'] = all_stat_df_caskey['Time_2'].astype(str) + '_' + all_stat_df_caskey['Participant']

#Label the viral load data
all_stat_df_caskey['VL_1'] = all_stat_df_caskey['Time_Par1'].map(vl_dict)
all_stat_df_caskey['VL_2'] = all_stat_df_caskey['Time_Par2'].map(vl_dict)
all_stat_df_caskey['Ave_VL'] = all_stat_df_caskey[['VL_1', 'VL_2']].mean(axis=1)

print('The mean viral load is:')
print(round(np.mean(all_stat_df_caskey['Ave_VL'])))
print('The standard deviation of the viral load is:')
print(round(np.std(all_stat_df_caskey['Ave_VL'])))
