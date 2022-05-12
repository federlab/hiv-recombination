import sys
# sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
#for running on desktop
from matplotlib.ticker import MaxNLocator
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


#In this file I am going to figure out how many snps per dist x time are available in the Zanini data

#for running on desktop
dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini/analysis/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/zanini/05-05-2022/'

#where to find things we need
lociDir = dataDir + 'segregating_Loci/'
timepointFile = dataDir + 'DaysToTimepoints.csv'

#make a dictionary to convert timepoint to days
timepointDF = pd.read_csv(timepointFile)
timepointDict = {}
for curr_par in timepointDF['Participant'].unique():
    #format dataframe properly
    par_timepoints = timepointDF[timepointDF['Participant'] == curr_par]
    par_timepoints = par_timepoints.drop('Participant', axis = 1)
    par_timepoints = par_timepoints.set_index(' Timepoint')
    #convert it to a dictionary
    par_timepoints = par_timepoints.to_dict()
    par_timepoints = par_timepoints[' Day']
    timepointDict[curr_par] = par_timepoints

#loop through the files and plot the numbers of segregating loci
all_seg = []

for curr_file in os.listdir(lociDir):
    #get the participant and timepoint
    curr_par = curr_file.split('_')[1]
    curr_frag = curr_file.split('_')[2]
    curr_frag = curr_frag.split('.')[0]

    #get the timepoint dictionary for the current participant
    curr_timepoints = timepointDict[curr_par]
    curr_seg = pd.read_pickle(lociDir + curr_file)
    curr_seg['timepoint'] = list(map(int, curr_seg['timepoint']))
    curr_seg['Day'] = curr_seg['timepoint'].apply(lambda x: curr_timepoints[x] if x in curr_timepoints.keys() else None )
    curr_seg['Participant'] = curr_par
    curr_seg['Fragment'] = curr_frag
    curr_seg['num_timepoints'] = len(curr_seg['timepoint'].unique())
    all_seg.append(curr_seg)

all_seg = pd.concat(all_seg, ignore_index= True)
print(all_seg)

#plot # snps by fragment
myplot = sns.FacetGrid(all_seg, col = 'Fragment')
myplot.map_dataframe(sns.histplot, x = 'Day', data = all_seg)
plt.xlabel("Day")
plt.ylabel("Number of Segregating Sites")
plt.tight_layout()
plt.savefig(outDir + "segregating_loci_fragment.jpg")
plt.close()

#plot # snps by participant
myplot = sns.FacetGrid(all_seg, col = 'Participant')
myplot.map_dataframe(sns.histplot, x = 'Day', data = all_seg)
plt.xlabel("Day")
plt.ylabel("Number of Segregating Sites")
plt.tight_layout()
plt.savefig(outDir + "segregating_loci_participant.jpg")
plt.close()

#plot the mean number of SNPS per timepoint
time_par_df = all_seg.groupby(['timepoint','Participant'])['timepoint'].count()
time_par_df = time_par_df.reset_index(name='snp_counts')
print("The quantiles of the SNP distribution per timepoint are")
my_quantiles = time_par_df['snp_counts'].quantile([0.25,0.5,0.75])
print(my_quantiles)
ax = sns.histplot(x = 'snp_counts', bins = 20, data = time_par_df)
ax.yaxis.set_major_locator(MaxNLocator(integer=True))
plt.ylabel("Count")
plt.xlabel("Number of Segregating Sites")
plt.tight_layout()
plt.savefig(outDir + "segregating_loci_timepoint.jpg")
plt.close()

#plot the mean number of SNPS per timepoint
frag_par_df = all_seg.groupby(['Fragment','Participant'])['timepoint'].count()
frag_par_df = frag_par_df.reset_index(name='snp_counts')
ax = sns.histplot(x = 'snp_counts', bins = 20, data = frag_par_df)
ax.yaxis.set_major_locator(MaxNLocator(integer=True))
plt.ylabel("Count")
plt.xlabel("Number of Segregating Sites")
plt.tight_layout()
plt.savefig(outDir + "segregating_loci_fragment_5_included.jpg")
plt.close()

frag_par_df = frag_par_df[frag_par_df['Fragment'] != 'F5']
print("The quantiles of the SNP distribution per fragment are")
my_quantiles = frag_par_df['snp_counts'].quantile([0.25,0.5,0.75])
print(my_quantiles)
ax = sns.histplot(x = 'snp_counts', bins = 20, data = frag_par_df)
ax.yaxis.set_major_locator(MaxNLocator(integer=True))
plt.ylabel("Count")
plt.xlabel("Number of Segregating Sites")
plt.tight_layout()
plt.savefig(outDir + "segregating_loci_fragment.jpg")
plt.close()