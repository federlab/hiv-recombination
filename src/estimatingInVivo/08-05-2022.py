import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import zaniniUtil as zu
from matplotlib import rcParams

#Today I want to plot read coverage over distance between loci
dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/zanini/08-05-2022/'
coverage_df = []

#Loop through all of the directories to get the filtered genotypes
for curr_data in os.listdir(dataDir):
    #only get the data directories, not hidden files
    if curr_data[0] == '.':
        continue
    
    #get the genotypes
    filteredGenotypes = pd.read_pickle(dataDir + curr_data + "/analysis/FilteredGenotypes")

    #summarize the number of reads per pair of loci
    read_counts = filteredGenotypes.groupby(['Locus_1', 'Locus_2'])
    for name, group in read_counts:
        coverage_df.append([name[0], name[1], np.sum(group['Count'])])

coverage_df = pd.DataFrame(coverage_df, columns = ['Locus_1', 'Locus_2', 'Read_Count'])
coverage_df['Distance'] = coverage_df['Locus_2'] - coverage_df['Locus_1']
coverage_df['Log_Count'] = np.log10(coverage_df['Read_Count'])

dist_coverage = []
#summarize the average read count by distance
coverage_group = coverage_df.groupby('Distance')
for name, group in coverage_group:
    dist_coverage.append([name, np.mean(group['Read_Count'])])
dist_coverage = pd.DataFrame(dist_coverage, columns = ['Distance', 'Average_Count'])
print(coverage_df['Read_Count'])

#################### Plot the Distribution of read coverage ###################
sns.scatterplot(data = coverage_df, x = 'Locus_1', y = 'Locus_2', hue = 'Log_Count', alpha = 0.01)
plt.savefig(outDir + 'Read_Coverage.jpg')
plt.close()

sns.scatterplot(data = dist_coverage, x = 'Distance', y = 'Average_Count')
plt.savefig(outDir + 'Dist_Coverage.jpg')
plt.close()



