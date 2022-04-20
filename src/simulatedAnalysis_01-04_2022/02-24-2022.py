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
#Now that we can do analysis using the diagonal scanning method, we can plot
dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_02_08/mu1e-04_rho1e-05_Ne10000_M800_rep1/numpy/'
outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_02_08/mu1e-04_rho1e-05_Ne10000_M800_rep1/'

#for running on desktop
dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_02_08/mu1e-04_rho1e-05_Ne10000_M800_rep1/numpy/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_02_08/mu1e-04_rho1e-05_Ne10000_M800_rep1/'

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

    # calculate R^2 values for our locus pairs
    r2List, distList, supportList = r2.calculate_R2_pairCounts(
        coCounts_arr, segregatingLoci)
    curr_df = pd.DataFrame(r2List, columns=['r2'])
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
        ave_r2 = curr_window['r2'].mean()
        center = winStart + (WINSIZE/2)
        patient_aves.append([center, winStart, winEnd, ave_r2])

patient_aves = pd.DataFrame(patient_aves, columns = ['center','window_start', 'window_end', 'average'])

#plot the results for the current participant
sns.set(rc={'figure.figsize':(15,5)})
myplot = sns.scatterplot(x = 'dist', y = 'r2', hue = 'date', data = timepoints, alpha = 0.5)
myplot.legend(loc='center left', bbox_to_anchor=(1.25, 0.5), ncol=2)
sns.lineplot(x = 'center', y = 'average', data = patient_aves, linewidth = 3)
plt.ylim(-0.1,1.1)
plt.xlim(-10,max(timepoints['dist']))
plt.xlabel("Distance Between Loci")
plt.ylabel("R^2 Value")
plt.tight_layout()
plt.savefig(outDir + curr_time + "_window_" + str(WINSIZE))
plt.close()
