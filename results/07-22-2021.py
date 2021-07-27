import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
import os
import r2Analysis as r2
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import zaniniUtil as zu
from lmfit import minimize, Parameters, fit_report

######################### Helper Functions ####################################
def hill_weir(rho, n, x):
    """code adapted from https://eacooper400.github.io/gen8900/exercises/ld-2.html"""
    return (((10+(rho*x))/((2+(rho*x))*(11+(rho*x))))*(1+(((3+(rho*x))*(12+(12*(rho*x))+((rho*x)**2)))/(n*(2+(rho*x))*(11+(rho*x))))))

def residual(params, x, data):
    rho = params['rho']
    n = params['n']

    #this is the equation we are using for our fit
    model = hill_weir(rho, n, x)
    return data-model

def make_vl_df(loadDataDir):
    """Loop through all the patients and make a dataframe with the viral loads
    for each timepoint"""
    allLoads = []
    for currFile in os.listdir(loadDataDir):
        participant = (currFile.split('_')[1]).split('.')[0]
        currDF = pd.read_csv(loadDataDir + currFile, sep = '\t')
        currDF['participant'] = participant
        allLoads.append(currDF)
    return pd.concat(allLoads)

###############################################################################

#Today we are going to fit equations to our data to try and estimate rho 
dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/analysis/'
vlDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/'
outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/zanini/paper_filtering/fits/'
dataFiles = os.listdir(dataDir)

WINSIZE = 20
WINSTEP = 5

allStats = []

#get all of the viral load data
allLoads = make_vl_df(vlDir)
print(allLoads, file = sys.stderr)

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

    unique_dates = np.unique(currData['date'])

    #perform the fit separately for each timepoint
    for date in unique_dates:
        timepoint_df = currData[currData['date'] == date]

        #start by setting parameters for our fit
        params = Parameters()
        params.add('rho', value=0.001)
        params.add('n', value=len(timepoint_df.index), vary=False)

        out = minimize(residual, params, args=(timepoint_df['dist'],), kws={'data': timepoint_df['r_squared']})
        print(currFile.split('_')[1])
        print(fit_report(out))

        rho_estimate = out.params['rho']
        n_current = out.params['n']
        x_vals = list(range(0, max(timepoint_df['dist'])))
        fitted_vals = [hill_weir(rho_estimate, n_current, x) for x in x_vals]
        fit_df = pd.DataFrame(list(zip(x_vals, fitted_vals)), columns= ['x_vals', 'fitted_vals' ])

        #make a dataframe where we save the fit vs viral load for each timepoint.
        if n_current > 20:
            vl = allLoads[allLoads['participant'] == currFile.split('_')[1]]
            #check if the viral load was measured at this timepoint
            if date > len(vl.index):
                #if we don't have a viral load measurement, continue
                continue
            vl = vl.iloc[[int(date) - 1]]
            vl = list(vl['Viral load [virions/ml]'])[0]

            rho_and_vl = [vl, float(rho_estimate), int(n_current), date, currFile.split('_')[1], (currFile.split('_')[2]).split('.')[0]]
            allStats.append(rho_and_vl)

        #get the averages
        #now we can calculate our moving average
        windowStarts = range(0,int(max(timepoint_df['dist'])), WINSTEP)
        current_aves = []
        #subset the dataframe based on distance
        for i in windowStarts:
            winStart = i
            winEnd = i + WINSIZE
            #get all of the datapoints in our window
            curr_window = timepoint_df[timepoint_df['dist'].between(winStart, winEnd)]
            if not curr_window.empty:
                ave_r = curr_window['r_squared'].mean()
                center = winStart + (WINSIZE/2)
                current_aves.append([center, winStart, winEnd, ave_r])
        

        current_aves = pd.DataFrame(current_aves, columns = ['center','window_start', 'window_end', 'average'])

        #plot the results for the current participant and fragment
        sns.set(rc={'figure.figsize':(15,5)})
        myplot = sns.scatterplot(x = 'dist', y = 'r_squared', data = timepoint_df, alpha = 0.5)
        sns.lineplot(x = 'center', y = 'average', data = current_aves, linewidth = 3)
        sns.lineplot(x = 'x_vals', y = 'fitted_vals', data = fit_df, linewidth = 3)
        plt.ylim(-0.1,1.1)
        plt.xlim(-10,max(timepoint_df['dist']))
        plt.xlabel("Distance Between Loci")
        plt.ylabel("R^2 Value")
        plt.tight_layout()
        plt.savefig(outDir + currFile.split('_')[1] + "_window_" + str(WINSIZE) + (currFile.split('_')[2]).split('.')[0] + "_time_"+ str(date))
        plt.close()


#get the viral load and rho estimates
allStats = pd.DataFrame(allStats, columns = ['viral_load', 'rho', 'n_snps', 'date', 'participant', 'fragment'])
#plot rho vs viral load color by participant
myplot = sns.scatterplot(x = 'viral_load', y = 'rho', hue = 'participant', data = allStats, alpha = 0.5)
myplot.set(xscale="log")
plt.xlabel("Viral Load [Virions/ml]")
plt.ylabel("Estimated Rho")
plt.ylim(-0.1,1.1)
plt.tight_layout()
plt.savefig(outDir + "Rho_vs_ViralLoad_Participant")
plt.close()

#plot rho vs viral load color by fragment
myplot = sns.scatterplot(x = 'viral_load', y = 'rho', hue = 'fragment', data = allStats, alpha = 0.5)
myplot.set(xscale="log")
plt.xlabel("Viral Load [Virions/ml]")
plt.ylabel("Estimated Rho")
plt.ylim(-0.1,1.1)
plt.tight_layout()
plt.savefig(outDir + "Rho_vs_ViralLoad_Fragment")
plt.close()


