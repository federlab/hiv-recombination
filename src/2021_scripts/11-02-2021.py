import sys
#for cluster run
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
# #for running on desktop
# sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import neher
import os
import statsmodels.api as sm
from statsmodels.tools.tools import add_constant
from statsmodels.sandbox.regression.predstd import wls_prediction_std
from lmfit import minimize, Parameters, fit_report

#directories for cluster run
dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/analysis/neher/'
dayDir =  '/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/analysis/'
outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/zanini/neher_analysis/fits/'

# # for running on desktop
# dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini/analysis/neher/'
# dayDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini/analysis/'
# outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/zanini/neher_analysis/fits/'

# Today I am going to use Weighted least squares to fit curves to my neher and
# leitner analysis. In this file, I am going to fit all of the data at once.

######################### Helper Functions ####################################
def neher_leitner(c0, c1, c2, dDeltaT):
    """The function Neher and Leitner fit to determine recombination rate"""
    return (c0 + (c1 * (1 - np.exp(-c2 * dDeltaT))))

def residual(params, dDeltaT, data, rec_error):
    c0 = params['c0']
    c1 = params['c1']
    c2 = params['c2']

    #this is the equation we are using for our fit
    model = neher_leitner(c0, c1, c2, dDeltaT)

    resids = data - model
    weighted_resids = resids * (1 + rec_error)
    return weighted_resids
#################################################################################


#Next we need to set some things for the run
#We are only using fragments 1-4
fragment_list = ['F1','F2', 'F3', 'F4']
par_list = ['p1', 'p2','p3', 'p4', 'p5', 'p6', 'p7', 'p8', 'p9', 'p10', 'p11']

#values to truncate our fits at
trunc_vals = range(5000, 10000, 5000)

#distance bins
BINWIDTH = 250
MIN_BIN = 0
MAX_BIN = 60000
CUTOFF = 0.03
SUCCESS = 0.01
RUNNAME = str(CUTOFF) + '_' + str(SUCCESS) +  "fragments1-4_truncating_smaller_ranges_fewervals"

#create lists to store all of the results in
rec_dfs = []
mut_dfs = []

#Loop through all of the files and get their information.
for currfile in os.listdir(dataDir):
    #get the participant and the fragment
    curr_par = currfile.split('_')[1]
    curr_frag = currfile.split('_')[2]
    curr_frag = curr_frag.split('.')[0]

    #check if we want to analyze the data from this participant
    if curr_par not in par_list or curr_frag not in fragment_list:
        continue
    
    #read the file into our dataframe
    curr_df = pd.read_csv(dataDir + currfile, index_col= 0)
    curr_df['Participant'] = curr_par
    curr_df['Fragment'] = curr_frag

    #check which type of tests were conducted
    if currfile.split('_')[0] == 'Recombination':
        rec_dfs.append(curr_df)
    else: mut_dfs.append(curr_df)
mutation_df = pd.concat(mut_dfs)
recombination_df = pd.concat(rec_dfs)

#make a place to store all of the test results
all_frequencies_patients = []

#convert the units to basepairs x generation
recombination_df['Dist_x_Time'] = 0.5 * recombination_df['Dist_x_Time']
mutation_df['Dist_x_Time'] = 0.5 * mutation_df['Dist_x_Time']

#create our bins
for bin_start in range(MIN_BIN, MAX_BIN, BINWIDTH):
    bin_end = bin_start + BINWIDTH


    #make a place to store the frequencies of observing events
    mut_frequencies = []
    mut_successes = []
    recomb_frequencies = []
    recomb_successes = []
    mut_tests = []
    recomb_tests = []

    #get all of the datapoints in our bin
    curr_recomb = recombination_df[recombination_df['Dist_x_Time'].between(bin_start, bin_end)]
    curr_mut = mutation_df[mutation_df['Dist_x_Time'].between(bin_start, bin_end)]


    #Calculate the frequencies in each bin
    if curr_recomb.shape[0] > 0:
        recomb_true = curr_recomb[curr_recomb['Test_Passed'] == True]
        recomb_frequencies.append(recomb_true.shape[0]/curr_recomb.shape[0])
        recomb_successes.append(recomb_true.shape[0])
        recomb_tests.append(curr_recomb.shape[0])

    else: 
        recomb_frequencies.append(0)
        recomb_successes.append(0)
        recomb_tests.append(0)

    if curr_mut.shape[0] > 0:
        mut_true = curr_mut[curr_mut['Test_Passed'] == True]
        mut_frequencies.append(mut_true.shape[0]/curr_mut.shape[0])
        mut_successes.append(mut_true.shape[0])
        mut_tests.append(curr_mut.shape[0])
    else: 
        mut_frequencies.append(0)
        mut_successes.append(0)
        mut_tests.append(0)

    all_frequencies = pd.DataFrame(list(zip(mut_frequencies,recomb_frequencies)),
                            columns = ['mut_frequencies', 'recomb_frequencies'])
    all_frequencies['window'] = bin_start
    all_frequencies['Mutation Tests'] = mut_tests
    all_frequencies['Recombination Tests'] = recomb_tests
    all_frequencies['Mutation Successes'] = mut_successes
    all_frequencies['Recombination Successes'] = recomb_successes
    all_frequencies_patients.append(all_frequencies)

#now we can make a dataframe with all our results to plot
all_frequencies_patients = pd.concat(all_frequencies_patients)
#get the error bars for each set of tests
all_frequencies_patients['Mut Error'] = 1 / np.sqrt(all_frequencies_patients['Mutation Tests'])
all_frequencies_patients['Recomb Error'] =  1/ np.sqrt(all_frequencies_patients['Recombination Tests'])

#We are going to try running the fits on truncated 
fit_df = []
estimate_df = []
for curr_trunc in trunc_vals:
    trunc_data = all_frequencies_patients[all_frequencies_patients['window'] < curr_trunc]
    ########################### Fitting ###########################################
    #start by setting parameters for our fit
    params = Parameters()
    params.add('c0', min = 0, max = 1, value = 0.1)
    params.add('c1', min = 0, max = 1, value = 0.26, vary = False)
    params.add('c2', min = 0, max = 0.1, value = 0.00005, vary = False)
    # params.add('c0', min = 0, max = 1, value = 0.1)
    # params.add('c1', min = 0, max = 1, value = 0.5)
    # params.add('c2', min = 0, max = 0.1, value = 0.5)

    out = minimize(residual, params, args=(all_frequencies_patients['window'], ), kws={'data': trunc_data['recomb_frequencies'],
                                                                                        'rec_error' : trunc_data['Recomb Error']},
                                                                                        max_nfev = 200000)

    c0_estimate = out.params['c0']
    c1_estimate = out.params['c1']
    c2_estimate = out.params['c2']
    x_vals = list(range(0, max(trunc_data['window'])))
    fitted_vals = [neher_leitner(c0_estimate, c1_estimate, c2_estimate, x) for x in x_vals]
    fitted_vals_paper = [neher_leitner(0.1, 0.26, .0000439, x) for x in x_vals]
    curr_fit_data = pd.DataFrame(list(zip(x_vals, fitted_vals, fitted_vals_paper)), columns= ['x_vals', 'fitted_vals', 'fitted_vals_paper' ])
    curr_fit_data['Trunc_Value'] = curr_trunc
    fit_df.append(curr_fit_data)
    estimate_df.append([c0_estimate, c1_estimate, c2_estimate, curr_trunc])

    print(fit_report(out), file = sys.stderr)

fit_df = pd.concat(fit_df)
estimate_df = pd.Dataframe(estimate_df)
########################### Plotting ###########################################

#plot the mutation tests
sns.set(rc={'figure.figsize':(20,5)})
plt.errorbar(x = all_frequencies_patients['window'], y = all_frequencies_patients['mut_frequencies'],
         yerr = all_frequencies_patients['Mut Error'], xerr = None, ls = 'none', ecolor = 'gray')
sns.lineplot(x = 'window', y = 'mut_frequencies', data = all_frequencies_patients, color = 'gray')    
plt.errorbar(x = all_frequencies_patients['window'], y = all_frequencies_patients['recomb_frequencies'],
         yerr = all_frequencies_patients['Recomb Error'], xerr = None, ls = 'none', ecolor = 'red')
sns.lineplot(x = 'window', y = 'recomb_frequencies', data = all_frequencies_patients, color = 'red')  
sns.lineplot(x = 'x_vals', y = 'fitted_vals_paper', data = fit_df, color = 'cyan', label = 'Neher Fit')
for curr_trunc in trunc_vals:
    sns.lineplot(x = 'x_vals', y = 'fitted_vals', data = fit_df[fit_df['Trunc_Value'] == curr_trunc], color = 'black', label = curr_trunc)
plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
plt.ylim(-0.1,0.6)
plt.xlabel("Distance x Time [BP X Generation]")
plt.ylabel("Frequency")
plt.tight_layout()
plt.savefig(outDir + "allTogether" + RUNNAME + ".jpg")
plt.close()

# #plot the estimates
# sns.set(rc={'figure.figsize':(20,5)})
# sns.scatterplot(x = 'Trunc_Value', y = 'C0 Estimate', data = estimate_df)
# sns.scatterplot(x = 'Trunc_Value', y = 'C1 Estimate', data = estimate_df)
# sns.scatterplot(x = 'Trunc_Value', y = 'C2 Estimate', data = estimate_df)
# plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
# plt.ylim(-0.1,0.6)
# plt.xlabel("Distance x Time [BP X Generation]")
# plt.ylabel("Frequency")
# plt.tight_layout()
# plt.savefig(outDir + "coefficient_estimates" + RUNNAME + ".jpg")
# plt.close()
