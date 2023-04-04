import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import os
import seaborn as sns
import matplotlib.pyplot as plt
import zaniniUtil as zu
from matplotlib import rcParams

#################################################### Figure 2 #######################################################
#Today I am quickly plotting the viral loads to add to my poster

dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/'
outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/croiAnalysis/'

dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/croiAnalysis/'

available_files = os.listdir(dataDir)

#data frame of viral loads for each patient
viralLoadData = zu.make_viral_load_df(dataDir)

#now plot the viral loads
rcParams.update({'font.size': 18, 'figure.figsize':(12, 6)})
myplot = sns.lineplot(x = 'Days from infection', y = 'Viral load [virions/ml]', hue = 'Participant', data = viralLoadData,
                     errorbar = None)
myplot.legend(ncol = 6)
myplot.set(yscale = 'log')
plt.xlabel("Days from Infection")
plt.ylabel("Viral Load [virions/ml]")
plt.xlim(-10, 3000)
plt.tight_layout()
plt.savefig(outDir + "viralLoads", dpi = 300)
plt.close()

quit()