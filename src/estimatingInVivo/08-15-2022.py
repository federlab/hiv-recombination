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
from matplotlib import rcParams

THRESHOLD = 50000

#Today I am making a plot of the timepoints alondside their viral loads. This 
#will serve to illustrate how the data is grouped

#For running on desktop
dataDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/"
vlDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/"
outDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/zanini/08-15-2022/"

#Make the dataframe containg D' ratios
stat_df = zu.combine_drats(dataDir)
stat_df = zu.label_vl_drats(stat_df, vlDir)

stat_df = stat_df[stat_df['Fragment'] != 'F5']
stat_df = stat_df[stat_df['Time_Diff'].gt(50)]

#Only get values for which the viral load doesn't leave the boundaries of the 
#initial and final measurments
stat_df = stat_df[stat_df['Monotonic'] == True]

#make a facet grid with 11 plots
fig, axs = plt.subplots(3, 4, figsize=(10, 10), sharex= True, sharey= True)
fig.delaxes(axs[-1][-1])

curr_par = 1



#loop through each facet and plot the data for one participant
for i in range(3):
    for j in range(4):
        #we only need to do 11 plots instead of 12
        if i == 2 and j == 3:
            break

        if i == 2 and j == 0:
            axs[i,j].set_xlabel('Time (days)')
            axs[i,j].set_ylabel('Viral Load (copies/ml)')

        #get the data for the participant in this facet
        curr_par_data = stat_df[stat_df['Participant'] == "p" + str(curr_par)]
        # print("********** Participant " + str(curr_par) + " **********")
        # print(curr_par_data['Day_1'].unique())
        # print(curr_par_data['Day_2'].unique())

        #now loop through pairs of timepoints and plot them on the grid
        grouped_par = curr_par_data.groupby(by = ['Day_1', 'Day_2'])
        for name, group in grouped_par:

            time_1 = name[0]
            time_2 = name[1]

            vl_1 = group['VL_1'].unique()[0]
            vl_2 = group['VL_2'].unique()[0]

            ave_vl = group['Ave_VL'].unique()[0]
            ave_time = np.mean([time_1, time_2])
            if  ave_vl > THRESHOLD:
                my_color = "tab:orange"
                my_color2 = "peru"
            else: 
                my_color = "tab:blue"
                my_color2 = "navy"

            # axs[i,j].axline((time_1, vl_1), (time_2, vl_2), color = my_color)
            axs[i,j].plot([time_1, time_2], [vl_1, vl_2], color = my_color)
            axs[i,j].axhline(THRESHOLD, color = "gray", linestyle = "--")
            axs[i,j].set_title("p" + str(curr_par))

        for name, group in grouped_par:
            time_1 = name[0]
            time_2 = name[1]

            vl_1 = group['VL_1'].unique()[0]
            vl_2 = group['VL_2'].unique()[0]

            ave_vl = group['Ave_VL'].unique()[0]
            ave_time = np.mean([time_1, time_2])
            if  ave_vl > THRESHOLD:
                my_color = "tab:orange"
                my_color2 = "saddlebrown"
            else: 
                my_color = "tab:blue"
                my_color2 = "navy"
            
            axs[i,j].plot(ave_time, ave_vl, color = my_color2, marker = "D", markersize = 2)

        #update the participant
        curr_par += 1

# plt.yscale('log')
plt.savefig(outDir + 'new_vl_prototype.jpg')
plt.close()

#plot all of the points on top of each other

#now loop through pairs of timepoints and plot them on the grid
grouped_par = stat_df.groupby(by = ['Day_1', 'Day_2'])
for name, group in grouped_par:

    time_1 = name[0]
    time_2 = name[1]

    vl_1 = group['VL_1'].unique()[0]
    vl_2 = group['VL_2'].unique()[0]

    ave_vl = group['Ave_VL'].unique()[0]
    ave_time = np.mean([time_1, time_2])
    if  ave_vl > THRESHOLD:
        my_color = "tab:orange"
        my_color2 = "peru"
    else: 
        my_color = "tab:blue"
        my_color2 = "navy"

    plt.plot([time_1, time_2], [vl_1, vl_2], color = my_color)
    plt.axhline(THRESHOLD, color = "gray", linestyle = "--")

for name, group in grouped_par:

    time_1 = name[0]
    time_2 = name[1]

    vl_1 = group['VL_1'].unique()[0]
    vl_2 = group['VL_2'].unique()[0]

    ave_vl = group['Ave_VL'].unique()[0]
    ave_time = np.mean([time_1, time_2])
    if  ave_vl > THRESHOLD:
        my_color = "tab:orange"
        my_color2 = "saddlebrown"
    else: 
        my_color = "tab:blue"
        my_color2 = "navy"


    plt.plot(ave_time, ave_vl, color = my_color2, marker = "D", markersize = 2)

plt.xlabel('Time (days)')
plt.ylabel('Viral Load (copies/ml)')
plt.savefig(outDir + 'all_vls_prototype.jpg')
plt.close()