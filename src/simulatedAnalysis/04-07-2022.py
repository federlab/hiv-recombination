import sys
# sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
#For running on Desktop
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import os
import numpy as np
import pandas as pd
import seaborn as sns
import plot_neher as plne
import matplotlib.pyplot as plt
import zaniniUtil as zu

#I am going to use this script to make figures for the Harris Lab Meeting

#For running on cluster
dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_03_17/'
outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_03_17/'

#For running on desktop
dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_02_24/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_02_24/'

# if not os.path.isdir(outDir):
#     os.mkdir(outDir)

# ###############################################################################
# #Make dataframes to store the resulting frequencies in
# all_dt_freqs = []
# all_d_freqs = []
# estimate_df = []

# for curr_dataset in os.listdir(dataDir):
#     print(curr_dataset)
#     if curr_dataset.split('_')[0] != 'mu1e-05':
#         continue
#     curr_rho = curr_dataset.split('_')[1]
#     curr_rho = curr_rho.split('o')[1]

#     #Get the current data
#     recombination_df = pd.read_pickle(
#                         dataDir + curr_dataset + "/neher_res/recombination")
#     recombination_df['dist'] = recombination_df['Locus_2'] - \
#                                 recombination_df['Locus_1']
#     recombination_df['Dist_x_Time'] = (recombination_df['Curr_Timepoint'] - \
#                 recombination_df['Last_Timepoint']) * recombination_df['dist']

#     mutation_df = pd.read_pickle(
#                         dataDir + curr_dataset + "/neher_res/mutation")
#     mutation_df['Dist_x_Time'] = mutation_df['Dist_x_Time'].abs()

#     ###### Plot by Dist_x_Time ######
#     #make the set of bins
#     dtbinset = [(x, x + 2500) for x in range(0, 55000, 2500)]
#     dt_freqs = plne.bin_curve(
#         recombination_df, mutation_df, dtbinset, 'Dist_x_Time')
#     dt_freqs['Dataset'] = curr_dataset
#     dt_freqs['Rho'] = curr_rho
#     all_dt_freqs.append(dt_freqs)
#     print(dtbinset[0])
#     zero_mut = dt_freqs[dt_freqs['window'] == dtbinset[0][1]]
#     zero_mut = zero_mut['mut_frequencies'].tolist()[0]
#     print(zero_mut)
  

#     ###### Plot by dist ######    
#     #make the set of bins
#     dbinset = [(x, x + 50) for x in range(0, 500, 50)]
#     d_freqs = plne.bin_curve(
#         recombination_df, mutation_df, dbinset, 'dist')
#     d_freqs['Dataset'] = curr_dataset
#     d_freqs['Rho'] = curr_rho
#     all_d_freqs.append(d_freqs)

#     # ########################### Fitting #######################################
#     #perform our fitting
#     try:
#         # coeffs, fit_df = plne.run_neher_fit(c0_fixed = False, lower_bounds = [0,0,0],
#         #                             upper_bounds = [1,1,1], 
#         #                             initial = [0.13, 0.1804, 0.000114],
#         #                             test_results = dt_freqs)
#         coeffs, fit_df = plne.run_neher_fit(c0_fixed = zero_mut, lower_bounds = [0,0],
#                                     upper_bounds = [1,1], 
#                                     initial = [0.1804, 0.000114],
#                                     test_results = dt_freqs)
#     except ValueError as my_err:
#         print(my_err, file = sys.stderr)
#         coeffs = [0,0,0]
#         fit_df = pd.DataFrame(columns= ['x_vals', 'fitted_vals'], index = range(1))

#     curr_estimate = plne.estimate_recombination_rate(c0 = coeffs[0], c1 = coeffs[1], c2 = coeffs[2])
#     estimate_df.append([curr_estimate, curr_dataset, curr_rho])
#     print(estimate_df)

#     #Plot our frequencies with fits
#     sns.set(rc={'figure.figsize':(20,5)}, font_scale = 2) 
#     plt.errorbar(x = dt_freqs['window'], y = dt_freqs['recomb_frequencies'],
#             yerr = dt_freqs['Recomb Error'], xerr = None, ls = 'none', ecolor = 'red')
#     sns.lineplot(x = 'window', y = 'recomb_frequencies', data = dt_freqs, color = 'red', label = 'Recombination Tests', linestyle = '--')  
#     sns.lineplot(x = 'x_vals', y = 'fitted_vals', data = fit_df, color = 'black', label = 'Our Fit',linewidth = 2)
#     sns.lineplot(x = 'window', y = 'mut_frequencies', data = dt_freqs, color = 'blue', )
#     plt.errorbar(x = dt_freqs['window'], y = dt_freqs['mut_frequencies'],
#         yerr = dt_freqs['Mut Error'], xerr = None, ls = 'none', ecolor = 'blue')
#     plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
#     plt.ylim(-0.1,0.6)
#     plt.xlabel("Distance x Time [BP X Generation]")
#     plt.ylabel("Frequency")
#     plt.tight_layout()
#     plt.savefig(outDir + curr_dataset + "/mutation_intercept_neherResults.jpg")
#     plt.close()


# #Concatenate all of the frequencies together
# all_d_freqs = pd.concat(all_d_freqs, ignore_index = True )
# all_dt_freqs = pd.concat(all_dt_freqs, ignore_index= True)

estimate_df = [[0.4044080732307525, 'mu1e-05_rho0.001_Ne10000_M800_rep3', '0.001'], [6.394795379963177e-05, 'mu1e-05_rho2e-05_Ne10000_M800_rep10', '2e-05'], [3.306754157217316e-05, 'mu1e-05_rho1e-05_Ne10000_M800_rep8', '1e-05'], [0.0005288469812606526, 'mu1e-05_rho2e-04_Ne10000_M800_rep8', '2e-04'], [0.0002588956746829431, 'mu1e-05_rho1e-04_Ne10000_M800_rep4', '1e-04'], [np.nan, 'mu1e-05_rho2e-06_Ne10000_M800_rep9', '2e-06'], [4.114286450906866e-05, 'mu1e-05_rho1e-05_Ne10000_M800_rep10', '1e-05'], [6.834663380677528e-05, 'mu1e-05_rho2e-05_Ne10000_M800_rep7', '2e-05'], [5.362708300849257e-05, 'mu1e-05_rho1e-05_Ne10000_M800_rep4', '1e-05'], [3.126523841991085e-05, 'mu1e-05_rho2e-06_Ne10000_M800_rep2', '2e-06'], [2.133767616376019e-05, 'mu1e-05_rho2e-06_Ne10000_M800_rep5', '2e-06'], [0.0005242415225672068, 'mu1e-05_rho2e-04_Ne10000_M800_rep4', '2e-04'], [7.754549566708137e-05, 'mu1e-05_rho2e-05_Ne10000_M800_rep6', '2e-05'], [4.2098649815313946e-05, 'mu1e-05_rho1e-05_Ne10000_M800_rep3', '1e-05'], [2.4923495570374208e-05, 'mu1e-05_rho1e-05_Ne10000_M800_rep7', '1e-05'], [0.012149799288308775, 'mu1e-05_rho0.001_Ne10000_M800_rep8', '0.001'], [7.127400432356403e-05, 'mu1e-05_rho2e-05_Ne10000_M800_rep1', '2e-05'], [0.0002520002187887563, 'mu1e-05_rho1e-04_Ne10000_M800_rep10', '1e-04'], [0.08264525192893288, 'mu1e-05_rho0.001_Ne10000_M800_rep10', '0.001'], [7.75778399583739e-05, 'mu1e-05_rho2e-05_Ne10000_M800_rep3', '2e-05'], [6.276523977500086e-07, 'mu1e-05_rho2e-06_Ne10000_M800_rep4', '2e-06'], [np.nan, 'mu1e-05_rho2e-06_Ne10000_M800_rep7', '2e-06'], [4.0757686516191425e-05, 'mu1e-05_rho1e-05_Ne10000_M800_rep5', '1e-05'], [0.0005076008550414735, 'mu1e-05_rho2e-04_Ne10000_M800_rep7', '2e-04'], [0.00027903649815096154, 'mu1e-05_rho1e-04_Ne10000_M800_rep6', '1e-04'], [np.nan, 'mu1e-05_rho2e-06_Ne10000_M800_rep8', '2e-06'], [6.119576529494383e-05, 'mu1e-05_rho2e-05_Ne10000_M800_rep5', '2e-05'], [6.928974743894301e-08, 'mu1e-05_rho2e-06_Ne10000_M800_rep6', '2e-06'], [0.0266090460253163, 'mu1e-05_rho0.001_Ne10000_M800_rep2', '0.001'], [0.4131125171360321, 'mu1e-05_rho0.001_Ne10000_M800_rep7', '0.001'], [6.105359811804407e-05, 'mu1e-05_rho2e-05_Ne10000_M800_rep8', '2e-05'], [0.0002027386028359209, 'mu1e-05_rho1e-04_Ne10000_M800_rep2', '1e-04'], [3.315989140636631e-05, 'mu1e-05_rho1e-05_Ne10000_M800_rep9', '1e-05'], [6.225102078515421e-05, 'mu1e-05_rho2e-05_Ne10000_M800_rep2', '2e-05'], [0.0004942938748801653, 'mu1e-05_rho2e-04_Ne10000_M800_rep1', '2e-04'], [0.00036395521963918684, 'mu1e-05_rho2e-04_Ne10000_M800_rep2', '2e-04']]

estimate_df = pd.DataFrame(estimate_df, columns=["Est_Rho", 'Dataset', 'Sim_Rho'])

############################# Plotting Estimate Accuracy ######################
# #Plot our estimates against each other 
#make the rho values ints rather than strings
rhoDict = {"0.001" : 0.001,
            "1e-04" : 0.0001,
            "2e-04" : 0.0002,
            "1e-05" : 0.00001,
            "2e-05" : 0.00002,
            "2e-06" : 0.000002}

intRhoList = []
for entry in estimate_df['Sim_Rho']:
    intRhoList.append(rhoDict[entry])
estimate_df['Sim_int_rho'] = intRhoList

#add jitter to the rho values
x_vals = np.linspace(0.00001, 0.000105, 10)
def jitter(values,j):
    return values + np.random.normal(j,0.1,values.shape)
estimate_df['Sim_int_rho']

print(estimate_df)
sns.set(rc={'figure.figsize':(10,5)}, font_scale = 1)
print(estimate_df['Sim_int_rho'])


ax = sns.stripplot(x = 'Sim_Rho', y = 'Est_Rho', data = estimate_df, jitter = True,
    order = ["0.001", "2e-04", "1e-04", "2e-05", "1e-05", "2e-06"])

# distance across the "X" or "Y" stipplot column to span, in this case 40%
label_width = 0.4

for tick, text in zip(ax.get_xticks(), ax.get_xticklabels()):
    sample_name = text.get_text()  # "X" or "Y"

    # calculate the median value for all replicates of either X or Y
    rho_val = rhoDict[sample_name]

    # plot horizontal lines across the column, centered on the tick
    ax.plot([tick-label_width/2, tick+label_width/2], [rho_val, rho_val],
            lw=2, color='k')

plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
plt.xlabel("Simulation Value of Rho")
plt.ylabel("Estimated Value of Rho")
plt.ylim(0.0000001, 0.1)
plt.yscale('log')
plt.tight_layout()
plt.savefig(outDir + "mutation_intercept_comparedEstimates_stripplot_all_rhos.jpg")
plt.close()


