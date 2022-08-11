import sys
# #for running on cluster
# sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
#for running on desktop
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import os
import numpy as np
import pandas as pd
from statistics import median
import seaborn as sns
import matplotlib.pyplot as plt

#With this script I just want to get the average estimated number of templates
dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini/estimated_depths'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/zanini/02-08-2022/'


f1 = open(dataDir + '/templateNums.txt', 'w')
all_summaries = []

for currfile in os.listdir(dataDir):
    print(currfile)
    if currfile.split('.')[1] != "tsv":
        continue

    depth_df = pd.read_csv(dataDir +  "/" +currfile, sep = '\t')
    allTemplateNums = depth_df['templates approx'].tolist()
    f1 = depth_df['F1_coverage'].tolist()
    f2 = depth_df['F2_coverage'].tolist()
    f3 = depth_df['F3_coverage'].tolist()
    f4 = depth_df['F4_coverage'].tolist()
    f5 = depth_df['F5_coverage'].tolist()
    f6 = depth_df['F6_coverage'].tolist()

    curr_summary = pd.DataFrame(list(zip(allTemplateNums, f1, f2, f3, f4, f5, f6)), columns = ['templates approx', 'F1_coverage', 'F2_coverage', 'F3_coverage', 'F4_coverage', 'F5_coverage', 'F6_coverage'])
    all_summaries.append(curr_summary)

all_summaries = pd.concat(all_summaries)

print(np.median(all_summaries['templates approx']))
print(np.mean(all_summaries['templates approx']))

sns.set(rc={'figure.figsize':(10, 15)}, font_scale = 2)
fig, ax = plt.subplots(6, sharex=True, sharey=True)
[sns.scatterplot(data = all_summaries, x = 'templates approx', y = 'F' + str(x) + '_coverage', ax = ax[x-1] )\
    for x in range(1, 7)]
[ax[i].set_yscale('log') for i in range(6)]
[ax[i].set_xscale('log') for i in range(6)]
plt.savefig(outDir + 'templates_vs_coverage.jpg')
plt.close()