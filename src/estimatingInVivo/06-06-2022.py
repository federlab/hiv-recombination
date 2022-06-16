import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

dataDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/shankarrapa/"
outDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/shankarrapa/"

#This file plots the viral loads seen in the shankarrapa data
vl_data = pd.read_csv(dataDir + "viralLoads.tsv", sep = '\t')
vl_data = vl_data[vl_data['LogRNA'] != '.']
vl_data = vl_data[vl_data['LogRNA'] != 'neg']
print(list(map(float, vl_data['LogRNA'])))
vl_data['LogRNA'] = list(map(float, vl_data['LogRNA']))
ax = sns.lineplot(x = 'T to SC (yr.)', y = 'LogRNA', data = vl_data, hue = 'Patient')
ax.axhline(np.log10(50000), color = 'k')
plt.savefig(outDir + "Viral_Loads")
print(vl_data)