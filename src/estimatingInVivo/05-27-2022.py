import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import zaniniUtil as zu
import os

#Today I am quickly plotting the viral loads to add to my poster

dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/'
outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/zanini/poster/'

dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/zanini/poster_peqg/'

available_files = os.listdir(dataDir)

#data frame of viral loads for each patient
viralLoadData = zu.make_viral_load_df(dataDir)

#now plot the viral loads
sns.set(rc={'figure.figsize':(10,5)})
plt.rcParams.update({'font.size': 30})
myplot = sns.lineplot(x = 'Days from infection', y = 'Viral load [virions/ml]', hue = 'Participant', data = viralLoadData,
                     ci = None)
myplot.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
myplot.set(yscale = 'log')
plt.xlabel("Days from Infection")
plt.ylabel("Viral Load [virions/ml]")
plt.xlim(-10, 3000)
plt.tight_layout()
plt.savefig(outDir + "viralLoads", dpi = 200)
plt.close()
