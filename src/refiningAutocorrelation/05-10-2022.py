import sys
# sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
#for running on desktop
from matplotlib.ticker import MaxNLocator
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


#I'm going to calculate the median number of tests above the 0.2 threshold in
#the Zanini data

dataDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/"
outDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/zanini/05-10-2022/"

#loop through the files and plot the numbers of segregating loci
all_seg = []

#loop through all of the directories
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
        all_seg.append(curr_seg)
    

all_seg = pd.concat(all_seg, ignore_index= True)
print(all_seg) 

# #plot # snps by fragment
# myplot = sns.FacetGrid(all_seg, col = 'Fragment', row = )
# myplot.map_dataframe(sns.histplot, x = 'Day', data = all_seg)
# plt.xlabel("Day")
# plt.ylabel("Number of Segregating Sites")
# plt.tight_layout()
# plt.savefig(outDir + "segregating_loci_fragment.jpg")
# plt.close()

# #plot # snps by participant
# myplot = sns.FacetGrid(all_seg, col = 'Participant')
# myplot.map_dataframe(sns.histplot, x = 'Day', data = all_seg)
# plt.xlabel("Day")
# plt.ylabel("Number of Segregating Sites")
# plt.tight_layout()
# plt.savefig(outDir + "segregating_loci_participant.jpg")
# plt.close()

# #plot the mean number of SNPS per timepoint
# time_par_df = all_seg.groupby(['timepoint','Participant'])['timepoint'].count()
# time_par_df = time_par_df.reset_index(name='snp_counts')
# print("The quantiles of the SNP distribution per timepoint are")
# my_quantiles = time_par_df['snp_counts'].quantile([0.25,0.5,0.75])
# print(my_quantiles)
# ax = sns.histplot(x = 'snp_counts', bins = 20, data = time_par_df)
# ax.yaxis.set_major_locator(MaxNLocator(integer=True))
# plt.ylabel("Count")
# plt.xlabel("Number of Segregating Sites")
# plt.tight_layout()
# plt.savefig(outDir + "segregating_loci_timepoint.jpg")
# plt.close()

#plot the mean number of SNPS per timepoint
frag_par_df = all_seg.groupby(['Fragment','Participant']).size()
frag_par_df = frag_par_df.reset_index(name='snp_counts')
print(frag_par_df)
ax = sns.histplot(x = 'snp_counts', bins = 50, data = frag_par_df)
ax.yaxis.set_major_locator(MaxNLocator(integer=True))
plt.ylabel("Count")
plt.xlabel("Number of D Ratios")
plt.tight_layout()
plt.savefig(outDir + "num_d_ratio_tests.jpg")
plt.close()

# frag_par_df = frag_par_df[frag_par_df['Fragment'] != 'F5']
# print("The quantiles of the SNP distribution per fragment are")
my_quantiles = frag_par_df['snp_counts'].quantile([0.25,0.5,0.75])
print(my_quantiles)
# ax = sns.histplot(x = 'snp_counts', bins = 20, data = frag_par_df)
# ax.yaxis.set_major_locator(MaxNLocator(integer=True))
# plt.ylabel("Count")
# plt.xlabel("Number of Segregating Sites")
# plt.tight_layout()
# plt.savefig(outDir + "segregating_loci_fragment.jpg")
# plt.close()