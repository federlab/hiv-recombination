import sys
# sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
#for running on desktop
from matplotlib.ticker import MaxNLocator
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import os
import csv
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


#I'm going to calculate the median number of tests by viral load from
#the Zanini data

dataDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/"
vlDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/"
outDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/zanini/05-17-2022/"

#loop through the files and plot the numbers of segregating loci
all_d_rats = []

#loop through all of the directories with linkage files
for curr_dir in os.listdir(dataDir):
    if curr_dir[0] == '.':
        continue
    for curr_file in os.listdir(dataDir + curr_dir):
        #get the participant and timepoint
        curr_par = curr_dir.split('_')[0]
        curr_frag = curr_dir.split('_')[1]

        #get the timepoint dictionary for the current participant
        curr_seg = pd.read_pickle(dataDir + curr_dir + "/linkage/d_ratio")
        curr_seg['Participant'] = curr_par
        curr_seg['Fragment'] = curr_frag
        all_d_rats.append(curr_seg)
    

all_d_rats = pd.concat(all_d_rats, ignore_index= True)
all_d_rats['Day_1'] = 2 * all_d_rats['Time_1']
all_d_rats['Day_2'] = 2 * all_d_rats['Time_2']
all_d_rats = all_d_rats[all_d_rats['Participant'] != 'p7']
all_d_rats = all_d_rats[all_d_rats['Participant'] != 'p4']


labeled_rats = []
#label each participants ratios with the corresponding viral loads
for curr_file in os.listdir(vlDir):
    dict_from_csv = {}
    if curr_file[0] == '.':
        continue

    #make a dictionary of timepoints and their viral loads
    with open(vlDir + curr_file, mode='r') as inp:
        reader = csv.reader(inp, delimiter= '\t')
        next(reader)
        dict_from_csv = {float(rows[0]):float(rows[1]) for rows in reader}
    
    #get the ratios for the participant
    participant = curr_file.split('.')[0]
    participant = participant.split('_')[1]
    curr_d_rats = all_d_rats[all_d_rats['Participant'] == participant]

    #label the ratios
    curr_d_rats = curr_d_rats[curr_d_rats['Day_1'].isin(dict_from_csv.keys())]
    curr_d_rats = curr_d_rats[curr_d_rats['Day_2'].isin(dict_from_csv.keys())]

    curr_d_rats['VL_1'] = curr_d_rats['Day_1'].map(lambda x: dict_from_csv[x])
    curr_d_rats['VL_2'] = curr_d_rats['Day_2'].map(lambda x: dict_from_csv[x])
    labeled_rats.append(curr_d_rats)

labeled_rats = pd.concat(labeled_rats, ignore_index= True)
labeled_rats['Ave_VL'] = labeled_rats[['VL_1', 'VL_2']].mean(axis=1)
print(labeled_rats)
del all_d_rats

#plot the mean number of SNPS per timepoint
ax = sns.histplot(x = 'Ave_VL', bins = 50, data = labeled_rats)
ax.yaxis.set_major_locator(MaxNLocator(integer=True))
# plt.yscale("log")
plt.ylabel("Number of D Ratios")
plt.xlabel("Average Viral Load")
plt.tight_layout()
plt.savefig(outDir + "num_d_ratio_tests_no_log.jpg")
plt.close()


my_quantiles = labeled_rats['Ave_VL'].quantile([0.25,0.5,0.75])
print(my_quantiles)