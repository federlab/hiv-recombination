import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
import os
import r2Analysis as r2
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import zaniniUtil as zu

#Today we are going to reconstitute the data we saved to make
#plots that are filtered matching what Zanini et al. did.
dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/analysis/'
outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/zanini/paper_filtering/'
dataFiles = os.listdir(dataDir)
WINSIZE = 20
WINSTEP = 5

allStats = []

for currFile in dataFiles:
    currData = pd.read_csv(dataDir + currFile)
    # our columns should be
    # Locus_1,Locus_2,p_A,p_B,AB_obs,Ab_obs,aB_obs,ab_obs,r_squared,d_prime,date
    #filter to get frequencies between 80 and 20 like Zanini did.
    currData = currData[currData['p_A'] < .8]
    currData = currData[currData['p_A'] > .2]
    currData = currData[currData['p_B'] < .8]
    currData = currData[currData['p_B'] > .2]
    currData['support'] = currData['AB_obs'] + currData['Ab_obs'] + currData['aB_obs'] + currData['ab_obs']
    #filter to get support of at least 200 reads
    currData = currData[currData['support'] > 200]
    #label the patient and fragment
    currData['participant'] = currFile.split('_')[1]
    currData['fragment'] = (currFile.split('_')[2]).split('.')[0]
    currData['dist'] = abs(currData['Locus_1'] - currData['Locus_2'])

    #get the averages
    #now we can calculate our moving average
    windowStarts = range(0,int(max(currData['dist'])), WINSTEP)
    current_aves = []
    #subset the dataframe based on distance
    for i in windowStarts:
        winStart = i
        winEnd = i + WINSIZE
        #get all of the datapoints in our window
        curr_window = currData[currData['dist'].between(winStart, winEnd)]
        if not curr_window.empty:
            ave_D = curr_window['d_prime'].mean()
            center = winStart + (WINSIZE/2)
            current_aves.append([center, winStart, winEnd, ave_D])
    
    current_aves = pd.DataFrame(current_aves, columns = ['center','window_start', 'window_end', 'average'])

    #plot the results for the current participant and fragment
    sns.set(rc={'figure.figsize':(15,5)})
    myplot = sns.scatterplot(x = 'dist', y = 'd_prime', hue = 'date', data = currData, alpha = 0.5)
    myplot.legend(loc='center left', bbox_to_anchor=(1.25, 0.5), ncol=2)
    sns.lineplot(x = 'center', y = 'average', data = current_aves, linewidth = 3)
    plt.ylim(-0.1,1.1)
    plt.xlim(-10,max(currData['dist']))
    plt.xlabel("Distance Between Loci")
    plt.ylabel("D' Value")
    plt.tight_layout()
    plt.savefig(outDir + currFile.split('_')[1] + "_window_" + str(WINSIZE) + (currFile.split('_')[2]).split('.')[0])
    plt.close()

    allStats.append(currData)



allStats = pd.concat(allStats)

#make the plots
fragmentList = ['F1', 'F2', 'F3', 'F4', 'F5', 'F6']

for currFrag in fragmentList:
    # get just the fragment
    plotData = allStats[allStats['fragment'] == currFrag]
    print(plotData, file = sys.stderr)

    #get the averages
    #now we can calculate our moving average
    windowStarts = range(0,int(max(plotData['dist'])), WINSTEP)
    patient_aves = []
    #subset the dataframe based on distance
    for i in windowStarts:
        winStart = i
        winEnd = i + WINSIZE
        #get all of the datapoints in our window
        curr_window = plotData[plotData['dist'].between(winStart, winEnd)]
        if not curr_window.empty:
            ave_D = curr_window['d_prime'].mean()
            center = winStart + (WINSIZE/2)
            patient_aves.append([center, winStart, winEnd, ave_D])
    
    patient_aves = pd.DataFrame(patient_aves, columns = ['center','window_start', 'window_end', 'average'])

    #plot the results for all our participants
    sns.set(rc={'figure.figsize':(15,5)})
    myplot = sns.scatterplot(x = 'dist', y = 'd_prime', hue = 'participant', data = plotData, alpha = 0.5)
    myplot.legend(loc='center left', bbox_to_anchor=(1.25, 0.5), ncol=2)
    sns.lineplot(x = 'center', y = 'average', data = patient_aves, linewidth = 3)
    plt.ylim(-0.1,1.1)
    plt.xlim(-10, max(plotData['dist']))
    plt.xlabel("Distance Between Loci")
    plt.ylabel("D' Value")
    plt.tight_layout()
    plt.savefig(outDir + "zanini_filtering_" + str(WINSIZE) + currFrag)
    plt.close()