import sys
# sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
#for running on desktop
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import os
import numpy as np
import neher
import plot_neher as plne
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import zaniniUtil as zu

#This file is essentially the same as 03-16-2022, but I am trying to use the functionality from plot neher
CUTOFF = 0.03

#For running on cluster
dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_02_24/mu1e-05_rho1e-04_Ne10000_M800_rep1'
outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_02_24/mu1e-05_rho1e-04_Ne10000_M800_rep1'

#For running on desktop
dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_02_24/mu1e-05_rho1e-04_Ne10000_M800_rep1'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_02_24/mu1e-05_rho1e-04_Ne10000_M800_rep1'

if not os.path.isdir(outDir):
    os.mkdir(outDir)


dataset_dist = []
dataset_dist_time = [] 

all_mut_df = []
all_rec_df = []

for curr_rep in range(1, 11): 

    print(curr_rep)
    recombination_df = pd.read_pickle(dataDir + "/neher_res/recombination")
    mutation_df = pd.read_pickle(dataDir + "/neher_res/mutation")

    currFile_gen = dataDir + "/analysis/FilteredGenotypes"
    currFile_loc = dataDir + "/analysis/FilteredLoci"

    filtered_genotypes = pd.read_pickle(currFile_gen)
    filtered_loci = pd.read_pickle(currFile_loc)

    #need to make timepoints strings for compatibility with neher analysis file
    filtered_genotypes['timepoint'] = list(map(str, filtered_genotypes['timepoint']))
    filtered_loci['timepoint'] = list(map(str, filtered_loci['timepoint']))
    mutation_df = neher.mutation_analysis(filtered_loci, filtered_loci['frag_len'].unique()[0], CUTOFF)
    mutation_df['Dist_x_Time'] = (mutation_df['Curr_Timepoint'] - mutation_df['Last_Timepoint']) * mutation_df['dist']
    all_mut_df.append(mutation_df)
    all_rec_df.append(recombination_df)

    #####################Distance x Time Plot##################  

    #the bins for our plot
    bin_option2 = [(x, x + 5000) for x in range(0, 100000, 2500)]
    all_frequencies_dt = plne.bin_curve(recombination_df, mutation_df, bin_option2, bin_type = 'Dist_x_Time')
    all_frequencies_dt['rep'] = curr_rep

    #####################Distance Plot##################

    #create our bins
    bin_set = [(0, 50), (50, 100), (100, 150), (150, 200), (200, 250), (250, 300), (300, 350), (350, 400), (400, 450), (450, 500) ]
    all_frequencies_d = plne.bin_curve(recombination_df, mutation_df, bin_set, bin_type = 'dist')
    all_frequencies_d['rep'] = curr_rep

    dataset_dist.append(all_frequencies_d)
    dataset_dist_time.append(all_frequencies_dt)

    break

all_mut_df = pd.concat(all_mut_df, ignore_index= True)
all_rec_df = pd.concat(all_rec_df, ignore_index= True)
dataset_dist = pd.concat(dataset_dist, ignore_index= True)
dataset_dist_time = pd.concat(dataset_dist_time, ignore_index= True)

sns.scatterplot(x = 'dist', y = 'Dist_x_Time', data = all_mut_df, alpha = 0.15, hue = 'Test_Passed')
plt.xlabel('Distance')
plt.ylabel("Distance x Time [BP X Generation]")
plt.savefig(outDir + "mutation_sanity")
plt.close()


sns.lineplot(x = 'window', y = 'mut_frequencies', data = dataset_dist_time, hue = 'rep')    
sns.lineplot(x = 'window', y = 'recomb_frequencies', data = dataset_dist_time, hue = 'rep', linestyle = '--' )
plt.ylim(0,0.6)
plt.xlabel("Distance x Time [BP X Generation]")
plt.ylabel("Frequency")
plt.tight_layout()
plt.savefig(outDir + "dist_time")
plt.close()


sns.lineplot(x = 'window', y = 'mut_frequencies', data = dataset_dist, hue = 'rep')    
plt.ylim(0,0.6)
plt.xlabel("Distance")
plt.ylabel("Frequency")
plt.tight_layout()
plt.savefig(outDir + "dist")
plt.close()