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


#For running on cluster
dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_02_24/'
outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_02_24/'

#For running on desktop
dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_02_24/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_02_24/'

if not os.path.isdir(outDir):
    os.mkdir(outDir)

# #make the set of bins

bin_option1 = [(0, 500), (500, 750), (750, 1000), (1000, 2000), (2000, 3000), (3000, 4000), 
    (4000, 5000), (5000, 7500), (7500, 10000), (10000, 15000), (15000, 20000), (20000, 30000), (30000, 40000), (40000, 50000) ,(50000, 60000)]
#bins that neher and leitner approximately used
bin_option2 = [(0,5000), (5000, 12500), (12500, 22500), (22500, 30000), (30000, 37500), (37500, 45000), (45000, 52500)]

bin_options = [[(x, x + 2500) for x in range(0, 55000, 2500)], bin_option1, bin_option2]
bin_names = ['2500 Bins', 'Finer Bins', '7500 Bins']

all_dt_freqs = []

for i in range(len(bin_options)):
    curr_bin = bin_options[i]
    curr_name = bin_names[i]

    ###############################################################################
    #Make dataframes to store the resulting frequencies in

    estimate_df = []

    for curr_dataset in os.listdir(dataDir):
        print(curr_dataset)
        if curr_dataset.split('_')[0] != 'mu1e-05':
            continue
        curr_rho = curr_dataset.split('_')[1]
        curr_rho = curr_rho.split('o')[1]

        #Get the current data
        recombination_df = pd.read_pickle(
                            dataDir + curr_dataset + "/neher_res/recombination")
        recombination_df['dist'] = recombination_df['Locus_2'] - \
                                    recombination_df['Locus_1']
        recombination_df['Dist_x_Time'] = (recombination_df['Curr_Timepoint'].values - \
                    recombination_df['Last_Timepoint'].values) * recombination_df['dist'].values

        mutation_df = pd.read_pickle(
                            dataDir + curr_dataset + "/neher_res/mutation")
        mutation_df['Dist_x_Time'] = mutation_df['Dist_x_Time'].abs()
        print(recombination_df.columns)
        print(np.mean(recombination_df['Supporting_Reads'].values))
        print(recombination_df['Supporting_Reads'].values.shape)
        print(recombination_df['pA'].values.shape)
        
        mpApB = np.mean(np.multiply(recombination_df['pB'].values, np.multiply(recombination_df['Supporting_Reads'].values, recombination_df['pA'].values)))
        print(mpApB)

        ###### Plot by Dist_x_Time ######
        dt_freqs = plne.bin_curve(
            recombination_df, mutation_df, curr_bin, 'Dist_x_Time')
        dt_freqs['Dataset'] = curr_dataset
        dt_freqs['Rho'] = curr_rho
        dt_freqs['bin_name'] = curr_name
        all_dt_freqs.append(dt_freqs)


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
        print(estimate_df)

        # #Plot our frequencies with fits
        # sns.set(rc={'figure.figsize':(20,5)}, font_scale = 2) 
        # plt.errorbar(x = dt_freqs['window'], y = dt_freqs['recomb_frequencies'],
        #         yerr = dt_freqs['Recomb Error'], xerr = None, ls = 'none', ecolor = 'red')
        # sns.lineplot(x = 'window', y = 'recomb_frequencies', data = dt_freqs, color = 'red', label = 'Recombination Tests', linestyle = '--')  
        # sns.lineplot(x = 'x_vals', y = 'fitted_vals', data = fit_df, color = 'black', label = 'Our Fit',linewidth = 2)
        # sns.lineplot(x = 'window', y = 'mut_frequencies', data = dt_freqs, color = 'blue', )
        # plt.errorbar(x = dt_freqs['window'], y = dt_freqs['mut_frequencies'],
        #     yerr = dt_freqs['Mut Error'], xerr = None, ls = 'none', ecolor = 'blue')
        # plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
        # plt.ylim(-0.1,0.6)
        # plt.xlabel("Distance x Time [BP X Generation]")
        # plt.ylabel("Frequency")
        # plt.tight_layout()
        # plt.savefig(outDir + "/" + curr_name + "_mutation_intercept_neherResults.jpg")
        # plt.close()
        break
   
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
    plt.savefig(outDir + "/" + curr_name + "_mutation_intercept_comparedEstimates_stripplot_all_rhos.jpg")
    plt.close()


all_dt_freqs = pd.concat(all_dt_freqs, ignore_index=True)
intRhoList = []
for entry in all_dt_freqs['Rho']:
    intRhoList.append(rhoDict[entry])
all_dt_freqs['Sim_int_rho'] = intRhoList

#plot the results for the current participant
myplot = sns.FacetGrid(all_dt_freqs, row = 'bin_name')
myplot.map_dataframe(sns.lineplot, x = 'window', y = 'recomb_frequencies', data = all_dt_freqs, hue = 'Sim_int_rho')
plt.xlabel("Window")
plt.ylim(-0.1,1.1)
plt.ylabel("P(Test Success)")
plt.tight_layout()
plt.savefig(outDir + "faceted_curves_neher.jpg")
plt.close()
