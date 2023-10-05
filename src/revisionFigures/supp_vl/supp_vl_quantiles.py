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


#For running on desktop
dataDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/"
vlDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/"
outDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/paper/supp_vl/"

#plot the estimates to show how accurate they are
params = {'figure.figsize':(7, 5), 'axes.labelsize': 6,'axes.titlesize':6,  
          'legend.fontsize': 6, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
          'legend.title_fontsize': 6}
rcParams.update(params)

#Make the dataframe containg D' ratios
stat_df = zu.combine_drats(dataDir)
stat_df = zu.label_vl_drats(stat_df, vlDir)

stat_df = stat_df[stat_df['Fragment'] != 'F5']
stat_df = stat_df[stat_df['Time_Diff'].gt(50)]
stat_df = stat_df[stat_df['Participant'] != 'p10']

#Only get values for which the viral load doesn't leave the boundaries of the 
#initial and final measurements
stat_df = stat_df[stat_df['Monotonic'] == True]

#Get the quantiles of the viral load distribution
stat_df['VL_Quantile'], bins = pd.qcut(stat_df['Ave_VL'], q = [0, 0.33, 0.66, 1], retbins= True, labels = ["Lower", "Middle", "Upper"])
print(bins)


fig, axs = plt.subplots(4,3, sharex = True, sharey = True)
ax_nums = [(0,0), (0,1), (0,2), (1,0), (1,1), (1,2), (2,0), (2,1), (2,2), (3,0), (3,1)]
par_list = stat_df['Participant'].unique()
par_list = sorted(par_list, key = lambda x: int(x.split("p")[1]))

#Make a different plot for each participant
for i in range(len(par_list)):
    curr_par = par_list[i]
    curr_df = stat_df[stat_df['Participant'] == curr_par]

    curr_axs_nums = ax_nums[i]
    ax0 = curr_axs_nums[0]
    ax1 = curr_axs_nums[1]

    #now loop through pairs of timepoints and plot them on the grid
    grouped_par = curr_df.groupby(by = ['Day_1', 'Day_2'])
    for name, group in grouped_par:
        
        time_1 = name[0]
        time_2 = name[1]

        vl_1 = group['VL_1'].unique()[0]
        vl_2 = group['VL_2'].unique()[0]

        ave_vl = group['Ave_VL'].unique()[0]
        ave_time = np.mean([time_1, time_2])
        if  ave_vl > bins[1] and ave_vl < bins[2]:
            my_color = "tab:gray"
            my_color2 = "dimgray"
        elif ave_vl > bins[2]:
            my_color = "tab:orange"
            my_color2 = "peru"
        else: 
            my_color = "tab:blue"
            my_color2 = "navy"

        axs[ax0][ax1].plot([time_1, time_2], [vl_1, vl_2], color = my_color, alpha = 0.5)
        # axs[ax0][ax1].axhline(bin, color = "red", linestyle = "--")

        axs[ax0][ax1].plot(ave_time, ave_vl, color = my_color2, marker = "D", markersize = 1)
        curr_par_num = curr_par.split("p")[1]
        axs[ax0][ax1].set_title("Participant " + curr_par_num)

    #Now plot the viral load trajectory for the current individual
    par_vls = pd.read_csv(vlDir + 'viralLoad_' + curr_par + '.tsv', sep = '\t',
            header = 0, names = ['Time', 'VL'])
    axs[ax0][ax1].plot(par_vls['Time'], par_vls['VL'], color = "k", linestyle = "--")

axs[3][0].set_xlabel('Time (days since EDI)')
axs[2][1].set_xlabel('Time (days since EDI)')
axs[2][2].set_xlabel('Time (days since EDI)')

axs[3][0].set_ylabel('Viral Load\n (copies/mL)')
axs[2][0].set_ylabel('Viral Load\n (copies/mL)')
axs[1][0].set_ylabel('Viral Load\n (copies/mL)')
axs[0][0].set_ylabel('Viral Load\n (copies/mL)')

plt.tight_layout()
plt.savefig(outDir + 'supp_vl' + '.jpg', dpi = 300)

plt.close()