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

GROUP_THRESHOLD_LIST = [7500, 10000, 25000, 50000, 100000, 200000]
NUM_BOOTSTRAPS = 1000
EXAMPLE_THRESHOLD = 50000

# I am going to make a two panel histogram. The first panel will be a histogram
# of the viral loads associated with timepoint pairs. The second panel will be
# a histogram of the viral load differences associated with LDMs.

#For running on desktop
dataDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/"
vlDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/"
outDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/revision/supp_figs/supp_vlHist/"

#For running on cluster
# dataDir = "/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/"
# vlDir = "/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/"
# outDir = "/net/feder/vol1/home/evromero/2021_hiv-rec/results/revision/supp_figs/supp_vlHist/"

###############################################################################
##################### Setting the figure parameters ###########################
###############################################################################
#plot the estimates to show how accurate they are
params = {'figure.figsize':(7, 3), 'axes.labelsize': 6,'axes.titlesize':6,  
          'legend.fontsize': 6, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
          'legend.title_fontsize': 6}
rcParams.update(params)

fig, axs = plt.subplots(1, 2, sharex= True)
plt.subplots_adjust(wspace = 0.3)

linewidth = 1

###############################################################################

#Make the dataframe containing D' ratios
stat_df = zu.combine_drats(dataDir)
stat_df = zu.label_vl_drats(stat_df, vlDir)

stat_df = stat_df[stat_df['Fragment'] != 'F5']
stat_df = stat_df[stat_df['Time_Diff'].gt(50)]

#Only get values for which the viral load doesn't leave the boundaries of the 
#initial and final measurments
stat_df = stat_df[stat_df['Monotonic'] == True]

#Group the dataframe by the timepoint pairs
pair_VLs = []
for name, group in stat_df.groupby(['Time_1', 'Participant', 'Time_2']):
    curr_VL = np.log10(group['Ave_VL'].unique()[0])
    pair_VLs.append(curr_VL)

sns.histplot(pair_VLs, ax = axs[0], bins = 20, color = 'k')
axs[0].set_ylabel('Number of Timepoint Pairs')
axs[0].set_xlabel('Log10(Average Viral Load (cp/mL))')

#Now make the second panel which has just viral load vs LDMs
stat_df['Log_Ave_VL'] = np.log10(stat_df['Ave_VL'])
axs[1].set_ylabel('Number of LDMs')
axs[1].set_xlabel('Log10(Average Viral Load (cp/mL))')

sns.histplot(stat_df['Log_Ave_VL'], ax = axs[1], bins = 20, color = 'k')

plt.savefig(outDir + 'vlHist.jpg', dpi = 300)

###############################################################################
# Additional analysis to get the median viral load for the lowest tertile
QUANTILE_LIST = [0, 0.33, 0.66, 1]

#Now print the average viral load for the lowest viral load bin
lower_label = "Lower \n<" + str(int(np.quantile(stat_df['Ave_VL'], 0.33)))
middle_label = "Middle \n" + str(int(np.quantile(stat_df['Ave_VL'], 0.33))) + " - " + str(int(np.quantile(stat_df['Ave_VL'], 0.66))) 
upper_label = "Upper \n>" + str(int(np.quantile(stat_df['Ave_VL'], 0.66)))

#First we need to group the data by the viral load quantiles
stat_df['VL_Quantile'], bins = pd.qcut(stat_df['Ave_VL'], q = QUANTILE_LIST, retbins= True, labels = [lower_label, middle_label, upper_label])

lowest_quantile = stat_df[stat_df['VL_Quantile'] == lower_label]
print('The median viral load for the lowest quantile is: ')
print(np.mean(lowest_quantile['Ave_VL']))
print ('The standard deviation of the viral load for the lowest quantile is: ')
print(np.std(lowest_quantile['Ave_VL']))







