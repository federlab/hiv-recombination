import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
#for running on desktop
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plot_neher as plne
import zaniniUtil as zu
from scipy import stats
from matplotlib import rcParams
from matplotlib import gridspec

NUM_BOOTSTRAPS = 1000

#In this file I am making a plot of the estimates for just three groups of data points
#For running on desktop
dataDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/"
vlDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/"
outDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/zanini/04-12-2023/"

#Make the dataframe containing D' ratios
stat_df = zu.combine_drats(dataDir)
stat_df = zu.label_vl_drats(stat_df, vlDir)

stat_df = stat_df[stat_df['Fragment'] != 'F5']
stat_df = stat_df[stat_df['Time_Diff'].gt(50)]

#Only get values for which the viral load doesn't leave the boundaries of the 
#initial and final measurments
stat_df = stat_df[stat_df['Monotonic'] == True]
print(min(stat_df['d_i']))

#First plot a histogram of viral loads
sns.histplot(x ='Ave_VL', data=stat_df)
plt.savefig(outDir+ 'vl_hist.png')

#Get the quantiles of the viral load distribution
vl_quantiles = stat_df['Ave_VL'].quantile([0, 0.33, 0.66, 1])
quantile_dict = {0:"0- .33", 1:".33 - .66", 2:".66 - 1"}
print(vl_quantiles)

my_quantiles = [(vl_quantiles[0], vl_quantiles[0.33]), (vl_quantiles[0.33], vl_quantiles[0.66]), (vl_quantiles[0.66], vl_quantiles[1])]

#Make a place to store the estimates
all_par_ests = []
group_size_df = []

#Make three even groups of D' ratios
for i in range(len(my_quantiles)):
    curr_quantile = my_quantiles[i]
    lower_thresh = curr_quantile[0]
    upper_thresh = curr_quantile[1]
    curr_stat_df = stat_df[stat_df['Ave_VL'].between(lower_thresh, upper_thresh, inclusive = "right")]

    print(len(curr_stat_df))
    #Get the current estimate
    lower_fit, upper_fit, estimate_df = plne.bootstrap_rho(curr_stat_df,
                                                            NUM_BOOTSTRAPS)
    estimate_df['Group'] = quantile_dict[i]
    
    #Add the size of the group to the dictionary which labels the axis
    group_size_df.append([curr_quantile, i, curr_stat_df.shape[0]])
    
    #Add all of the current results
    all_par_ests.append(estimate_df )

all_par_ests = pd.concat(all_par_ests, ignore_index = True)

sns.boxplot(x ='Group', y ='Estimated_Rho', data = all_par_ests, palette = sns.color_palette("rocket"))
plt.ylim(0, 0.0003)
plt.xlabel("Viral Load Quantile")
plt.ylabel("Estimated Recombination Rate " + r'$\hat{\rho}$')
plt.tight_layout()
plt.savefig(outDir+ 'vl_quantile_boxplot.png')
    