import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
#for running on desktop
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sim_downsamp import make_downsampled_read_dist
from matplotlib import rcParams
from matplotlib import gridspec
from scipy.stats import binned_statistic

DIST_SAMPLES = 100000

#Downsampling parameters based on the Zanini et al. sequencing scheme
MIN_READ_LEN = 400
MAX_READ_LEN = 700
HALF_READ_LEN = 300

#This file will plot the empirical distributions of coverage in the simulated
#data and in the in vivo data

#Today I want to plot read coverage over distance between loci
dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/paper/supp_coverage/'
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

# dist_coverage = []
# #summarize the average read count by distance
# coverage_group = coverage_df.groupby('Distance')
# for name, group in coverage_group:
#     dist_coverage.append([name, np.mean(group['Read_Count'])])
# dist_coverage = pd.DataFrame(dist_coverage, columns = ['Distance', 'Average_Count'])

binned_cov, binedges, bin_nums = binned_statistic(
                coverage_df['Distance'].to_numpy(), 
                coverage_df['Read_Count'].to_numpy(), bins = 100)
invivo_cov = pd.DataFrame({'distance': binedges[:-1], 'coverage': binned_cov})

sim_coverage = make_downsampled_read_dist(DIST_SAMPLES, MIN_READ_LEN, MAX_READ_LEN, HALF_READ_LEN)
print(sim_coverage)
binned_cov, binedges, bin_nums = binned_statistic(
                sim_coverage.index.to_numpy(), 
                sim_coverage['coverage'].to_numpy(), bins = 100)
sim_cov = pd.DataFrame({'distance': binedges[:-1], 'coverage': binned_cov})

#################### Plot the Distribution of read coverage ###################
params = {'font.size': 8, 'axes.labelsize': 8,'axes.titlesize':8,  'legend.fontsize': 8, 'xtick.labelsize': 8, 'ytick.labelsize': 8}
plt.rcParams.update(params)
linewidth = 1

fig, axs = plt.subplots(1,2, figsize = (7, 3.5), sharex = True)

sns.lineplot(data = invivo_cov, x = 'distance', y = 'coverage', ax = axs[0], color = 'k')
axs[0].set_ylabel('Mean Number of Reads Overlapping Loci')
axs[0].set_xlabel(r'Distance Between Loci' + r' (bp)')

sns.lineplot(data = sim_cov, x = 'distance', y = 'coverage', ax = axs[1], color = 'k')
axs[1].set_ylabel('Proportion of Reads Overlapping Loci')
axs[1].set_xlabel(r'Distance Between Loci' + r' (bp)')

plt.savefig(outDir + 'Dist_Coverage.jpg', dpi = 300)
plt.close()