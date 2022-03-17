import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
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
dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_02_24/'
outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_02_24/'

#For running on desktop
dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_02_24/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_02_24/'

if not os.path.isdir(outDir):
    os.mkdir(outDir)

###############################################################################
#Make dataframes to store the resulting frequencies in
all_dt_freqs = []
all_d_freqs = []
estimate_df = []

i = 1

for curr_dataset in os.listdir(dataDir):
    curr_rho = curr_dataset.split('_')[1]
    curr_rho = curr_rho.split('o')[1]
    if not os.path.isdir(outDir + curr_dataset):
        os.mkdir(outDir + curr_dataset)

    #Get the current data
    recombination_df = pd.read_pickle(
                        dataDir + curr_dataset + "/neher_res/recombination")
    recombination_df['dist'] = recombination_df['Locus_2'] - \
                                recombination_df['Locus_1']
    recombination_df['Dist_x_Time'] = (recombination_df['Curr_Timepoint'] - \
                recombination_df['Last_Timepoint']) * recombination_df['dist']

    mutation_df = pd.read_pickle(
                        dataDir + curr_dataset + "/neher_res/mutation")
    mutation_df['Dist_x_Time'] = mutation_df['Dist_x_Time'].abs()

    ###### Plot by Dist_x_Time ######
    #make the set of bins
    dtbinset = [(x, x + 2500) for x in range(0, 55000, 2500)]
    dt_freqs = plne.bin_curve(
        recombination_df, mutation_df, dtbinset, 'Dist_x_Time')
    dt_freqs['Dataset'] = curr_dataset
    dt_freqs['Rho'] = curr_rho
    all_dt_freqs.append(dt_freqs)

    ###### Plot by dist ######    
    #make the set of bins
    dbinset = [(x, x + 50) for x in range(0, 500, 50)]
    d_freqs = plne.bin_curve(
        recombination_df, mutation_df, dbinset, 'dist')
    d_freqs['Dataset'] = curr_dataset
    d_freqs['Rho'] = curr_rho
    all_d_freqs.append(d_freqs)

    # ########################### Fitting #######################################
    #perform our fitting
    try:
        coeffs, fit_df = plne.run_neher_fit(c0_fixed = False, lower_bounds = [0,0,0],
                                    upper_bounds = [1,1,1], 
                                    initial = [0.13, 0.1804, 0.000114],
                                    test_results = dt_freqs)
    except ValueError as my_err:
        print(my_err, file = sys.stderr)
        coeffs = [0,0,0]
        fit_df = pd.DataFrame(columns= ['x_vals', 'fitted_vals'], index = range(1))

    curr_estimate = plne.estimate_recombination_rate(c0 = coeffs[0], c1 = coeffs[1], c2 = coeffs[2])
    estimate_df.append([curr_estimate, curr_dataset, curr_rho])

    # ########################### Plotting ###########################################
    #Plot our frequencies with fits
    sns.set(rc={'figure.figsize':(20,5)}, font_scale = 2)
    # plt.errorbar(x = dt_freqs['window'], y = dt_freqs['mut_frequencies'],
    #         yerr = dt_freqs['Mut Error'], xerr = None, ls = 'none', ecolor = 'gray')
    # sns.lineplot(x = 'window', y = 'mut_frequencies', data = dt_freqs, color = 'gray', label = 'Mutation Tests')    
    plt.errorbar(x = dt_freqs['window'], y = dt_freqs['recomb_frequencies'],
            yerr = dt_freqs['Recomb Error'], xerr = None, ls = 'none', ecolor = 'red')
    sns.lineplot(x = 'window', y = 'recomb_frequencies', data = dt_freqs, color = 'red', label = 'Recombination Tests')  
    # sns.lineplot(x = 'x_vals', y = 'fitted_vals_paper', data = fit_df, color = 'blue', label = 'Neher Fit', linewidth = 2, linestyle = '--')
    sns.lineplot(x = 'x_vals', y = 'fitted_vals', data = fit_df, color = 'black', label = 'Our Fit',linewidth = 2)
    plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
    plt.ylim(-0.1,0.6)
    plt.xlabel("Distance x Time [BP X Generation]")
    plt.ylabel("Frequency")
    plt.tight_layout()
    plt.savefig(outDir + curr_dataset + "/neherResults.jpg")
    plt.close()

    #Plot our frequencies without fits
    sns.set(rc={'figure.figsize':(20,5)}, font_scale = 2)
    plt.errorbar(x = dt_freqs['window'], y = dt_freqs['mut_frequencies'],
            yerr = dt_freqs['Mut Error'], xerr = None, ls = 'none', ecolor = 'gray')
    sns.lineplot(x = 'window', y = 'mut_frequencies', data = dt_freqs, color = 'gray', label = 'Mutation Tests')    
    plt.errorbar(x = dt_freqs['window'], y = dt_freqs['recomb_frequencies'],
            yerr = dt_freqs['Recomb Error'], xerr = None, ls = 'none', ecolor = 'red')
    sns.lineplot(x = 'window', y = 'recomb_frequencies', data = dt_freqs, color = 'red', label = 'Recombination Tests')  
    # sns.lineplot(x = 'x_vals', y = 'fitted_vals_paper', data = fit_df, color = 'blue', label = 'Neher Fit', linewidth = 2, linestyle = '--')
    # sns.lineplot(x = 'x_vals', y = 'fitted_vals', data = fit_df, color = 'black', label = 'Our Fit',linewidth = 2)
    plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
    plt.ylim(-0.1,0.6)
    plt.xlabel("Distance x Time [BP X Generation]")
    plt.ylabel("Frequency")
    plt.tight_layout()
    plt.savefig(outDir + curr_dataset + "/neherResults_nofit.jpg")
    plt.close()

    if i == 10:
        break
    print(i)
    i += 1