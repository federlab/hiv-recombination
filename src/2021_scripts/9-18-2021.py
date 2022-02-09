import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
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


available_files = os.listdir(dataDir)

#data frame of viral loads for each patient
viralLoadData = zu.make_viral_load_df(dataDir)

#now plot the viral loads
sns.set(rc={'figure.figsize':(10,5)})
myplot = sns.lineplot(x = 'Days from infection', y = 'Viral load [virions/ml]', hue = 'Participant', data = viralLoadData,
                     ci = None)
myplot.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
myplot.set(yscale = 'log')
plt.xlabel("Days from Infection")
plt.ylabel("Viral Load [virions/ml]")
plt.xlim(-10, 3000)
plt.tight_layout()
plt.savefig(outDir + "viralLoads")
plt.close()

#make a second plot but size dots by count
# sns.set(rc={'figure.figsize':(15,5)})
# myplot = sns.scatterplot(x = 'window', y = 'mut_frequencies', hue = 'Participant', data = all_frequencies_patients, alpha = 0.5,
#                             palette = 'Greys', size = 'Mutation Tests')
# myplot = sns.scatterplot(x = 'window', y = 'recomb_frequencies', hue = 'Participant', data = all_frequencies_patients, alpha = 0.5,
#                             size = 'Recombination Tests')
# myplot = sns.lineplot(x = 'window', y = 'recomb_frequencies', hue = 'Participant', data = all_frequencies_patients, alpha = 0.5,
#                     ci = None)
# myplot.legend(loc='center left', bbox_to_anchor=(1.25, 0.5), ncol=2)
# plt.ylim(-0.1,1.1)
# plt.xlim(-10, MAX_BIN)
# plt.xlabel("Distance Between Loci")
# plt.ylabel("Frequency")
# plt.tight_layout()
# plt.savefig(outDir + "neher_scatter" + RUNNAME )
# plt.close()