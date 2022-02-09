import sys
# #for running on cluster
# sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
#for running on desktop
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import os
import numpy as np
import neher
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl


#Today I am continuing work with the Slim Dataset.
#This script takes in the neher results and plots them
DATASET_LIST = ['r1', 'r2', 'r3']

estimate_df = []

# #For running on desktop
dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/slimDatasets/'

# #For running on cluster
# dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/slimDatasets/'
# outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/slimDatasets/'

enclosing_dir = '2022_01_25/'

all_datasets_freqs = []

# #Today I am going to be using the data that Alison simulated using SLIM
# #This file will take the text files and make them into numpy arrays that are
# #saved
# for curr_dataset in os.listdir(dataDir + enclosing_dir):
#     print(curr_dataset)
#     if curr_dataset.split('-')[0] != "mu1e":
#         continue
#     #make a directory for output 
#     if not os.path.isdir(outDir + enclosing_dir + curr_dataset):
#         os.mkdir(outDir + enclosing_dir + curr_dataset)
#     currDir = dataDir + enclosing_dir + curr_dataset + "/neher_res/"

#     #Start by getting the dataframes with the results of the mutation and recombination tests
#     mutation_df = pd.read_pickle(currDir + "mutation")
#     recombination_df = pd.read_pickle(currDir + "recombination")

#     #create our bins
#     bin_option1 = [(0,250), (250, 500), (500, 750), (750, 1000), (1000, 2000), (2000, 3000), (3000, 4000), 
#     (4000, 5000), (5000, 7500), (7500, 10000), (10000, 15000), (15000, 20000), (20000, 30000), (30000, 40000), (40000, 50000) ,(50000, 60000)]
#     #bins that neher and leitner approximately used
#     bin_option2 = [(0,5000), (5000, 12500), (12500, 22500), (22500, 30000), (30000, 37500), (37500, 45000), (45000, 52500)]
#     bin_lists = [bin_option1, bin_option2]

#     for i in range(len(bin_lists)):
#         custom_bins = bin_lists[i]
#         #Get the test results for the relevant bins
#         all_frequencies = []


#         for currBin in custom_bins:
#             bin_start = currBin[0]
#             bin_end = currBin[1]
#             #added 1/14/22 plot the bins at the center instead of the start
#             bin_center = int((currBin[1] - currBin[0])/2) + bin_start

#             #make a place to store the frequencies of observing events
#             mut_frequencies = []
#             mut_successes = []
#             recomb_frequencies = []
#             recomb_successes = []
#             mut_tests = []
#             recomb_tests = []

#             #get all of the datapoints in our bin
#             curr_recomb = recombination_df[recombination_df['Dist_x_Time'].between(bin_start, bin_end)]
#             curr_mut = mutation_df[mutation_df['Dist_x_Time'].between(bin_start, bin_end)]


#             #Calculate the frequencies in each bin
#             if curr_recomb.shape[0] > 0:
#                 recomb_true = curr_recomb[curr_recomb['Test_Passed'] == True]
#                 recomb_frequencies.append(recomb_true.shape[0]/curr_recomb.shape[0])
#                 recomb_successes.append(recomb_true.shape[0])
#                 recomb_tests.append(curr_recomb.shape[0])

#             else: 
#                 recomb_frequencies.append(0)
#                 recomb_successes.append(0)
#                 recomb_tests.append(0)

#             if curr_mut.shape[0] > 0:
#                 mut_true = curr_mut[curr_mut['Test_Passed'] == True]
#                 mut_frequencies.append(mut_true.shape[0]/curr_mut.shape[0])
#                 mut_successes.append(mut_true.shape[0])
#                 mut_tests.append(curr_mut.shape[0])
#             else: 
#                 mut_frequencies.append(0)
#                 mut_successes.append(0)
#                 mut_tests.append(0)

#             curr_frequencies = pd.DataFrame(list(zip(mut_frequencies,recomb_frequencies)),
#                                 columns = ['mut_frequencies', 'recomb_frequencies'])
#             curr_frequencies['window'] = bin_end
#             curr_frequencies['Mutation Tests'] = mut_tests
#             curr_frequencies['Recombination Tests'] = recomb_tests
#             curr_frequencies['Mutation Successes'] = mut_successes
#             curr_frequencies['Recombination Successes'] = recomb_successes
#             all_frequencies.append(curr_frequencies)

#         all_frequencies = pd.concat(all_frequencies)
#         all_frequencies['Sim_Rho'] = curr_dataset.split("_")[1]
#         all_frequencies['bin'] = i
#         all_frequencies['rep'] = curr_dataset.split("_")[-1]
#         all_frequencies = all_frequencies.reset_index()
#         all_datasets_freqs.append(all_frequencies)

#         # ########################### Fitting #######################################
#         #get the error bars for each set of tests
#         all_frequencies['Mut Error'] = 1 / np.sqrt(all_frequencies['Mutation Tests'])
#         all_frequencies['Recomb Error'] =  1/ np.sqrt(all_frequencies['Recombination Tests'])

#         #perform our fitting
#         coeffs, fit_df = neher.run_neher_fit(c0_fixed = False, lower_bounds = [0,0,0],
#                                     upper_bounds = [1,1,1], 
#                                     initial = [0.13, 0.1804, 0.000114],
#                                     test_results = all_frequencies)

#         curr_estimate = neher.estimate_recombination_rate(c0 = coeffs[0], c1 = coeffs[1], c2 = coeffs[2])
#         estimate_df.append([curr_estimate, i, curr_dataset.split("_")[1]])
#         # ########################### Plotting ###########################################
#         #Plot our frequencies with fits
#         sns.set(rc={'figure.figsize':(20,5)}, font_scale = 2)
#         plt.errorbar(x = all_frequencies['window'], y = all_frequencies['mut_frequencies'],
#                 yerr = all_frequencies['Mut Error'], xerr = None, ls = 'none', ecolor = 'gray')
#         sns.lineplot(x = 'window', y = 'mut_frequencies', data = all_frequencies, color = 'gray', label = 'Mutation Tests')    
#         plt.errorbar(x = all_frequencies['window'], y = all_frequencies['recomb_frequencies'],
#                 yerr = all_frequencies['Recomb Error'], xerr = None, ls = 'none', ecolor = 'red')
#         sns.lineplot(x = 'window', y = 'recomb_frequencies', data = all_frequencies, color = 'red', label = 'Recombination Tests')  
#         # sns.lineplot(x = 'x_vals', y = 'fitted_vals_paper', data = fit_df, color = 'blue', label = 'Neher Fit', linewidth = 2, linestyle = '--')
#         sns.lineplot(x = 'x_vals', y = 'fitted_vals', data = fit_df, color = 'black', label = 'Our Fit',linewidth = 2)
#         plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
#         plt.ylim(-0.1,0.6)
#         plt.xlabel("Distance x Time [BP X Generation]")
#         plt.ylabel("Frequency")
#         plt.tight_layout()
#         plt.savefig(outDir + enclosing_dir + curr_dataset + "/neherResults_" + str(i) +  ".jpg")
#         plt.close()

#         #Plot our frequencies with fits
#         sns.set(rc={'figure.figsize':(20,5)}, font_scale = 2)
#         plt.errorbar(x = all_frequencies['window'], y = all_frequencies['mut_frequencies'],
#                 yerr = all_frequencies['Mut Error'], xerr = None, ls = 'none', ecolor = 'gray')
#         sns.lineplot(x = 'window', y = 'mut_frequencies', data = all_frequencies, color = 'gray', label = 'Mutation Tests')    
#         plt.errorbar(x = all_frequencies['window'], y = all_frequencies['recomb_frequencies'],
#                 yerr = all_frequencies['Recomb Error'], xerr = None, ls = 'none', ecolor = 'red')
#         sns.lineplot(x = 'window', y = 'recomb_frequencies', data = all_frequencies, color = 'red', label = 'Recombination Tests')  
#         # sns.lineplot(x = 'x_vals', y = 'fitted_vals_paper', data = fit_df, color = 'blue', label = 'Neher Fit', linewidth = 2, linestyle = '--')
#         # sns.lineplot(x = 'x_vals', y = 'fitted_vals', data = fit_df, color = 'black', label = 'Our Fit',linewidth = 2)
#         plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
#         plt.ylim(-0.1,0.6)
#         plt.xlabel("Distance x Time [BP X Generation]")
#         plt.ylabel("Frequency")
#         plt.tight_layout()
#         plt.savefig(outDir + enclosing_dir + curr_dataset + "/neherResults_nofit_" + str(i) +  ".jpg")
#         plt.close()

# estimate_df = pd.DataFrame(estimate_df, columns = ['Estimate', 'Binning', 'Sim_Rho'])
# all_datasets_freqs = pd.concat(all_datasets_freqs)

# estimate_df.to_csv(outDir + enclosing_dir + "/estimate_df")
# all_datasets_freqs.to_csv(outDir + enclosing_dir + "/all_dataset_freqs_df")

estimate_df = pd.read_csv(outDir + enclosing_dir + "/estimate_df")
all_datasets_freqs = pd.read_csv(outDir + enclosing_dir + "/all_dataset_freqs_df")

# ########################### Plotting ###########################################
#Plot our frequencies with fits
binset0 = all_datasets_freqs[all_datasets_freqs['bin'] == 0]
sns.set(rc={'figure.figsize':(20,5)}, font_scale = 1)
myplot = sns.FacetGrid(all_datasets_freqs, col = 'Sim_Rho', row = 'bin')
myplot.map_dataframe(sns.lineplot, x = 'window', y = 'mut_frequencies', hue = 'rep' , linestyle = '--')    
myplot.map_dataframe(sns.lineplot, x = 'window', y = 'recomb_frequencies', hue = 'rep')  
# plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
plt.ylim(-0.1,0.6)
plt.xlabel("Distance x Time [BP X Generation]")
plt.ylabel("Frequency")
plt.tight_layout()
plt.savefig(outDir + enclosing_dir + "/neherResults_binsTogether_" + str(0) +  ".jpg")
plt.close()

# #Plot our frequencies with fits
# binset1 = all_datasets_freqs[all_datasets_freqs['bin'] == 1]
# sns.set(rc={'figure.figsize':(20,5)}, font_scale = 1)
# myplot = sns.FacetGrid(binset1, col = 'Sim_Rho')
# myplot.map_dataframe(sns.lineplot, x = 'window', y = 'mut_frequencies', hue = 'rep', linestyle = '--')    
# myplot.map_dataframe(sns.lineplot, x = 'window', y = 'recomb_frequencies', hue = 'rep')  
# # plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
# plt.ylim(-0.1,0.6)
# plt.xlabel("Distance x Time [BP X Generation]")
# plt.ylabel("Frequency")
# plt.tight_layout()
# plt.savefig(outDir + enclosing_dir + "/neherResults_" + str(1) +  ".jpg")
# plt.close()

# #Plot our estimates against each other 
#make the rho values ints rather than strings
rhoDict = {"rho0" : 1e-20,
            "rho0.000105" : 0.000105,
            "rho1e-05" : 0.00001}
intRhoList = []
for entry in estimate_df['Sim_Rho']:
    intRhoList.append(rhoDict[entry])
estimate_df['Sim_int_rho'] = intRhoList
binning_dict = { 0: "Finer",
                1 : "Coarser"}
bin_name_list = []
for entry in estimate_df["Binning"]:
    bin_name_list.append(binning_dict[entry])
estimate_df['Bin Name'] = bin_name_list
#add jitter to the rho values
x_vals = np.linspace(0.00001, 0.000105, 10)
def jitter(values,j):
    return values + np.random.normal(j,0.1,values.shape)
estimate_df['Sim_int_rho']

print(estimate_df)
sns.set(rc={'figure.figsize':(10,5)}, font_scale = 1)

not_zeros = estimate_df[estimate_df['Sim_int_rho'] != 1e-20]
# myplot = sns.FacetGrid(estimate_df, col = 'Sim_int_rho', row = 'Binning', sharex = False, sharey = False)
# myplot.map_dataframe(sns.scatterplot, x = 'Sim_int_rho', y = 'Estimate', hue = 'rep')    
myplot = sns.stripplot(x = 'Sim_Rho', y = 'Estimate', data = not_zeros, hue = 'Bin Name', jitter = True)
sns.pointplot(x = 'Sim_Rho', y = 'Sim_int_rho', data = not_zeros, capsize = 0.7, ci = 0, scale = 0, color = 'black')
# plt.plot(x_vals, x_vals, '-r', label = 'y=x')
plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
myplot.set(yscale = "log")
plt.xlabel("Simulation Value of Rho")
plt.ylabel("Estimated Value of Rho")
plt.tight_layout()
plt.savefig(outDir + enclosing_dir + "/comparedEstimates_stripplot.jpg")
plt.close()

