import sys
# sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
#for running on desktop
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plot_neher as plne
from scipy import optimize

#In this script I am looking at the intercept and asymptotic behavior of the 
#autocorrelation curves as a function of rho

#I tried measuring the autocorrelation of the D statistic, but it looks like we
#are getting a lot of noise. So I am going to try setting up an an initial 
#thresholding value to only run tests after high linkage is initially seen.
THRESHOLD = 0.2

#For running on cluster
dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_02_24/'
outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_02_24/'

dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_02_24/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_02_24/autocorr_fine_tuning/'

estimate_df = [] 


for curr_data in os.listdir(dataDir):
    print(curr_data, file = sys.stderr)
    curr_data = "mu1e-05_rho1e-05_Ne10000_M800_rep8"

    #only get the data directories, not hidden files
    if curr_data[0] == '.':
        continue
    run_info = curr_data.split('_')
    sim_rho = run_info[1]
    sim_rho = sim_rho[3:]
    rep = run_info[-1]

    #make a place to store our output
    currOut = outDir + curr_data
    if not os.path.exists(currOut):
        os.mkdir(currOut)
    
    linkage_file = dataDir + curr_data + "/linkage/r2_and_D"

    #first just try and print the array
    rd_arr = pd.read_pickle(linkage_file)
    rd_arr.dropna(inplace = True)
    rd_arr['support'] = rd_arr['AB_obs'] + rd_arr['Ab_obs'] + rd_arr['aB_obs'] + rd_arr['ab_obs']
    rd_arr = rd_arr.loc[(rd_arr['AB_obs'] > 0) & (rd_arr['aB_obs'] > 0)  & (rd_arr['Ab_obs'] > 0)  & (rd_arr['ab_obs'] > 0) ]

    #group by the pair of loci
    grouped_loci = rd_arr.groupby(['Locus_1', 'Locus_2'])

    #make a dataframe to save the results in
    stat_df = []

    #loop over the pairs of loci to calculate -log(D_(t_i+1)/D_(t_i))
    for name, group in grouped_loci:
        #loop over the timepoints
        group_times = group['timepoint'].unique()
        group_times.sort()
        for i in range(len(group_times) - 1):
            curr_time = group_times[i]
            next_time = group_times[i+1]

            d_i = group[group['timepoint'] == curr_time]
            d_i = d_i['d_prime'].tolist()[0]
            d_i_1 = group[group['timepoint'] == next_time]
            d_i_1 = (d_i_1['d_prime']).tolist()[0]

            if d_i < THRESHOLD:
                # print('D_i is zero')
                continue
            curr_val = -np.log(d_i_1/d_i)
            stat_df.append([curr_val, name[0], name[1], next_time - curr_time, d_i])


    stat_df = pd.DataFrame(stat_df, columns = ['d_ratio', 'Locus_1', 'Locus_2', 'Time_Diff', 'd_i'])
    stat_df['Dist_X_Time'] = (stat_df['Locus_2'] - stat_df['Locus_1']) * stat_df['Time_Diff']
    stat_df = stat_df[stat_df['d_ratio'].between(-10,10)]
    stat_df = stat_df[stat_df['Dist_X_Time'].between(0, 50000)]
    stat_df.replace([np.inf, -np.inf], np.nan, inplace=True)
    stat_df.dropna(inplace = True)

    coeffs, fit_dat = optimize.curve_fit(plne.neher_leitner, stat_df['Dist_X_Time'], stat_df['d_ratio'], p0 = [0, 0.26, .0000439])
    x_vals = stat_df['Dist_X_Time'].unique()
    fit_vals = [plne.neher_leitner(x, coeffs[0], coeffs[1], coeffs[2]) for x in x_vals]
    sns.scatterplot(x = 'Dist_X_Time', y = 'd_ratio', data = stat_df, alpha = 0.05, hue = 'd_i')
    sns.lineplot(x = 'Dist_X_Time', y = 'd_ratio', data = stat_df, estimator = np.mean)
    sns.lineplot(x = x_vals, y = fit_vals)
    plt.savefig(currOut + "/auto_plot_50k.jpg")
    plt.close()
    print(coeffs)
    print(coeffs[1] * coeffs[2])
    estimate_df.append([coeffs[0], coeffs[1], coeffs[2], coeffs[1] * coeffs[2], curr_data, sim_rho])  

    coeffs, fit_dat = optimize.curve_fit(plne.neher_leitner_c0_fixed, stat_df['Dist_X_Time'], stat_df['d_ratio'], p0 = [0.26, .0000439])
    x_vals = stat_df['Dist_X_Time'].unique()
    fit_vals = [plne.neher_leitner_c0_fixed(x, coeffs[0], coeffs[1]) for x in x_vals]
    sns.scatterplot(x = 'Dist_X_Time', y = 'd_ratio', data = stat_df, alpha = 0.05, hue = 'd_i')
    sns.lineplot(x = 'Dist_X_Time', y = 'd_ratio', data = stat_df, estimator = np.mean)
    sns.lineplot(x = x_vals, y = fit_vals)
    plt.savefig(currOut + "/auto_plot_50k_fixed.jpg")
    plt.close()
    print(coeffs)
    print(coeffs[0] * coeffs[1])
    estimate_df.append([0, coeffs[0], coeffs[1], coeffs[0] * coeffs[1], curr_data, sim_rho])  
    break


estimate_df = pd.DataFrame(estimate_df, columns=["C0", "C1", "C2", "Est_Rho", 'Dataset', 'Sim_Rho'] )
estimate_df["C0_+_C1"] = estimate_df['C0'] + estimate_df['C1']
print(estimate_df, file = sys.stderr)
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

sns.set(rc={'figure.figsize':(10,5)}, font_scale = 1)
# print(estimate_df['Sim_int_rho'])


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
plt.ylim(0.000000001, 0.1)
plt.yscale('log')
plt.tight_layout()
plt.savefig(outDir + "autocorr_comparedEstimates_stripplot_50k.jpg")
plt.close()

#first plot

sns.set(rc={'figure.figsize':(10,5)}, font_scale = 1)


ax = sns.stripplot(x = 'Sim_Rho', y = 'C1', data = estimate_df, jitter = True,
    order = ["0.001", "2e-04", "1e-04", "2e-05", "1e-05", "2e-06"])
plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
plt.xlabel("Simulation Value of Rho")
plt.ylabel("Estimated Value of C1")
plt.tight_layout()
plt.savefig(outDir + "autocorr_C1_50k.jpg")
plt.close()

ax = sns.stripplot(x = 'Sim_Rho', y = 'C0', data = estimate_df, jitter = True,
    order = ["0.001", "2e-04", "1e-04", "2e-05", "1e-05", "2e-06"])
plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
plt.xlabel("Simulation Value of Rho")
plt.ylabel("Estimated Value of C0")
plt.tight_layout()
plt.savefig(outDir + "autocorr_C0_50k.jpg")
plt.close()

ax = sns.stripplot(x = 'Sim_Rho', y = 'C2', data = estimate_df, jitter = True,
    order = ["0.001", "2e-04", "1e-04", "2e-05", "1e-05", "2e-06"])
plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
plt.xlabel("Simulation Value of Rho")
plt.ylabel("Estimated Value of C2")
plt.tight_layout()
plt.savefig(outDir + "autocorr_C2_50k.jpg")
plt.close()