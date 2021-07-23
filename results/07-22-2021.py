import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
import os
import r2Analysis as r2
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import zaniniUtil as zu
from lmfit import minimize, Parameters

######################### Helper Functions ####################################

def residual(params, x, data):
    rho = params['rho']
    n = params['n']

    #code adapted from https://eacooper400.github.io/gen8900/exercises/ld-2.html
    model = (((10+(rho*x))/((2+(rho*x))*(11+(rho*x))))*(1+(((3+(rho*x))*(12+(12*(rho*x))+((rho*x)**2)))/(n*(2+(rho*x))*(11+(rho*x))))))
    return data-model

###############################################################################

#Today we are going to fit equations to our data to try and estimate rho 
dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/analysis/'
outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/zanini/paper_filtering/'
dataFiles = os.listdir(dataDir)

allStats = []

#We can start by making 
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

    #start by setting parameters for our fit
    params = Parameters()
    params.add('rho', value=0.001, vary=False)
    params.add('n', value=len(currData.index), vary=False)

    out = minimize(residual, params, args=(currData['dist'],), kws={'data': currData['r_squared']})
    out.params.pretty_print()

    break

#     #plot the results for the current participant and fragment
#     sns.set(rc={'figure.figsize':(15,5)})
#     myplot = sns.scatterplot(x = 'dist', y = 'd_prime', hue = 'date', data = currData, alpha = 0.5)
#     myplot.legend(loc='center left', bbox_to_anchor=(1.25, 0.5), ncol=2)
#     sns.lineplot(x = 'center', y = 'average', data = current_aves, linewidth = 3)
#     plt.ylim(-0.1,1.1)
#     plt.xlim(-10,max(currData['dist']))
#     plt.xlabel("Distance Between Loci")
#     plt.ylabel("D' Value")
#     plt.tight_layout()
#     plt.savefig(outDir + currFile.split('_')[1] + "_window_" + str(WINSIZE) + (currFile.split('_')[2]).split('.')[0])
#     plt.close()

#     allStats.append(currData)

# allStats = pd.concat(allStats)

#Now, we need to do the fit
#code from https://eacooper400.github.io/gen8900/exercises/ld-2.html
#for this small test we can subset the data to just one example

