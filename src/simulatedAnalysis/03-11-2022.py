import sys
# sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
#For running on Desktop
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import os
import r2Analysis as r2
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import zaniniUtil as zu

#Today I am debugging and checking that the slim simulations worked out

#for running on desktop
dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_02_24/mu1e-05_rho1e-04_Ne10000_M800_rep8/numpy/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_02_24/mu1e-05_rho1e-04_Ne10000_M800_rep8/'

#parameters for our moving median
#how much to advance the window start by
WINSTEP = 5
#the size of the window
WINSIZE = 20
WINSIZE_ALL = 10


time_list = ['t1', 't2', 't3', 't4', 't5', 't6', 't7', 't8']

timepoints = []

#loop through the timepoints
for curr_time in time_list:
    print(curr_time)
    #load our current counts
    coCounts_arr = np.load(dataDir + 'slim_formatted_' + curr_time + '.npy')

    #find the segregating sites
    segregatingLoci = zu.find_segregating_diagonal(coCounts_arr)

    #continue if there are no segregating sites
    if segregatingLoci.empty:
        continue

    # calculate D values for our locus pairs
    dList, distList, supportList = r2.calculate_R2_pairCounts(
        coCounts_arr, segregatingLoci, statistic = 'D')
    curr_df = pd.DataFrame(dList, columns=['D'])
    curr_df['dist'] = distList
    curr_df['support'] = supportList
    curr_df['date'] = str(curr_time)
    timepoints.append(curr_df)

#put the data from all the timepoints into a nice dataframe
timepoints = pd.concat(timepoints)
timepoints = timepoints.sort_values(by=['date'])

#now we can calculate our moving average
windowStarts = range(0,int(max(timepoints['dist'])), WINSTEP)
patient_aves = []

##############THRESHOLD = 50 supporting reads ###############
timepoints= timepoints[timepoints['support']>=50]

#subset the dataframe based on distance
for i in windowStarts:
    winStart = i
    winEnd = i + WINSIZE
    #get all of the datapoints in our window
    curr_window = timepoints[timepoints['dist'].between(winStart, winEnd)]
    if not curr_window.empty:
        ave_r2 = curr_window['D'].mean()
        center = winStart + (WINSIZE/2)
        patient_aves.append([center, winStart, winEnd, ave_r2])

patient_aves = pd.DataFrame(patient_aves, columns = ['center','window_start', 'window_end', 'average'])

#plot the results for the current participant
sns.set(rc={'figure.figsize':(15,5)})
myplot = sns.FacetGrid(timepoints, row = 'date')
myplot.map_dataframe(sns.scatterplot, x = 'dist', y = 'D', data = timepoints, alpha = 0.5)
plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5), ncol=2)
sns.lineplot(x = 'center', y = 'average', data = patient_aves, linewidth = 3)
plt.ylim(-0.1,1.1)
plt.xlim(-10,max(timepoints['dist']))
plt.xlabel("Distance Between Loci")
plt.ylabel("D Value")
plt.tight_layout()
plt.savefig(outDir + curr_time + "_window_" + str(WINSIZE))
plt.close()
