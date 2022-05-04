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
import autocorrelation as autocorr
from scipy import optimize
from scipy import stats

#In this script I am looking at the intercept and asymptotic behavior of the 
#autocorrelation curves as a function of rho

#I tried measuring the autocorrelation of the D statistic, but it looks like we
#are getting a lot of noise. So I am going to try setting up an an initial 
#thresholding value to only run tests after high linkage is initially seen.
THRESHOLD = 0.2
j = 0

#For running on cluster
dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_04_20/'
outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_04_20/'

dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_04_20/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_04_20/Autocorrelation/'

estimate_df = [] 

distTimeMax = [50000]

for curr_data in os.listdir(dataDir):
    print(curr_data)

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
    d_ratio_file = dataDir + curr_data + "/linkage/d_ratio"

    #04-25-2022 saving data for Alison
    # linkage_df = pd.read_pickle(linkage_file)
    # linkage_df.to_csv(outDir + curr_data + "/r2_and_D.csv")
    # break

    if os.path.exists(d_ratio_file):
        stat_df = pd.read_pickle(d_ratio_file)
    else:
        stat_df = autocorr.calculate_d_ratios(linkage_file, THRESHOLD)
        stat_df.to_pickle(d_ratio_file)

    # sns.scatterplot(x = 'Dist_X_Time', y = 'd_i', data = stat_df, alpha = 0.05)
    # sns.lineplot(x = 'Dist_X_Time', y = 'd_i', data = stat_df, estimator = np.median, color = "orange")
    # plt.savefig(currOut + "/d_i_median.jpg")



    for i in distTimeMax:
        stat_df_curr = stat_df[stat_df['Dist_X_Time'].between(0, i)]
        x_vals = stat_df_curr['Dist_X_Time'].unique()
        binned_rat, binedges, bin_nums = stats.binned_statistic(stat_df_curr['Dist_X_Time'].to_numpy(), stat_df_curr['d_ratio'].to_numpy(), bins = 100)
        print(binned_rat)
        print(binedges)


        #Fit the data and add our estimate to the dataframe
        # coeffs, fit_dat = optimize.curve_fit(plne.neher_leitner, stat_df_curr['Dist_X_Time'], stat_df_curr['d_ratio'], p0 = [0, 0.26, .0000439])
        # fit_vals = [plne.neher_leitner(x, coeffs[0], coeffs[1], coeffs[2]) for x in x_vals]
        sns.scatterplot(x = 'Dist_X_Time', y = 'd_ratio', data = stat_df_curr, alpha = 0.05, hue = 'd_i')
        sns.lineplot(x = 'Dist_X_Time', y = 'd_ratio', data = stat_df_curr, estimator = np.mean)
        sns.lineplot(x = binedges[:-1], y = binned_rat)
        # sns.lineplot(x = x_vals, y = fit_vals)
        plt.savefig(currOut + "/auto_plot_binned" + str(i) +".jpg")
        plt.close()

        coeffs, fit_dat = optimize.curve_fit(plne.neher_leitner, binedges[:-1], binned_rat, p0 = [0, 0.26, .0000439])
        fit_vals = [plne.neher_leitner(x, coeffs[0], coeffs[1], coeffs[2]) for x in x_vals]
        # sns.scatterplot(x = 'Dist_X_Time', y = 'd_ratio', data = stat_df_curr, alpha = 0.05, hue = 'd_i')
        # sns.lineplot(x = 'Dist_X_Time', y = 'd_ratio', data = stat_df_curr, estimator = np.mean)
        # sns.lineplot(x = binedges[:-1], y = binned_rat)
        # sns.lineplot(x = x_vals, y = fit_vals)
        # plt.savefig(currOut + "/auto_plot_" + str(i) +".jpg")
        # plt.close()

        estimate_df.append([coeffs[0], coeffs[1], coeffs[2], coeffs[1] * coeffs[2], curr_data, sim_rho, i, 'curve', curr_data])  
        print("final estimate is")
        print(estimate_df)



estimate_df = pd.DataFrame(estimate_df, columns=["C0", "C1", "C2", "Est_Rho", 'Dataset', 'Sim_Rho', 'distTimeMax', 'style', 'data'] )


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



# for curr_rho in np.unique(intRhoList):
#     sns.set(rc={'figure.figsize':(10,5)}, font_scale = 1)
#     curr_estimate_df = estimate_df[estimate_df['Sim_int_rho'] == curr_rho]

#     ax = sns.stripplot(x = 'style', y = 'Est_Rho', data = curr_estimate_df, jitter = True, hue = 'distTimeMax',
#         order = ["line", "curve"])
#     # distance across the "X" or "Y" stipplot column to span, in this case 40%
#     label_width = 0.4

#     for tick, text in zip(ax.get_xticks(), ax.get_xticklabels()):
#         sample_name = text.get_text()  # "X" or "Y"

#         # calculate the median value for all replicates of either X or Y
#         rho_val = curr_rho

#         # plot horizontal lines across the column, centered on the tick
#         ax.plot([tick-label_width/2, tick+label_width/2], [rho_val, rho_val],
#                 lw=2, color='k')
#     plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
#     plt.ylim(0.000000001, 0.1)
#     plt.yscale('log')
#     plt.xlabel("Fit Style")
#     plt.ylabel("Estimated Value of Rho")
#     plt.tight_layout()
#     plt.savefig(outDir + "autocorr_fit_style_rho_" + str(curr_rho) + ".jpg")
#     plt.close()

# # print(estimate_df['Sim_int_rho'])


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
plt.savefig(outDir + "autocorr_comparedEstimates_stripplot_50k_binned_fit.jpg")
plt.close()

# ax = sns.scatterplot(x = 'Sim_int_rho', y = 'Est_Rho', data = estimate_df)
# ax = sns.scatterplot(x = 'Sim_int_rho', y = 'Sim_int_rho', data = estimate_df, label = "y = x")
# plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
# plt.xlabel("Simulation Value of Rho")
# plt.ylabel("Estimated Value of Rho")
# # plt.ylim(0.000000001, 0.1)
# plt.yscale('log')
# plt.xscale('log')
# plt.tight_layout()
# plt.savefig(outDir + "autocorr_comparedEstimates_scatter_5k.jpg")
# plt.close()

# #fitting a line to the error
# estimate_df['error'] = np.log10(estimate_df['Est_Rho']) - np.log10(estimate_df['Sim_int_rho'])
# def line_func(x, m, b):
#     return m*x + b

# coeffs, covs = optimize.curve_fit(line_func, np.log10(estimate_df['Sim_int_rho']), estimate_df['error'])
# print("Line fit to error is")
# print(coeffs)
# x_vals = np.unique(intRhoList)
# fit_vals = [line_func(np.log10(x), coeffs[0], coeffs[1]) for x in x_vals]

# print(estimate_df, file = sys.stderr)
# ax = sns.scatterplot(x = 'Sim_int_rho', y = 'error', data = estimate_df)
# sns.lineplot(x = x_vals, y = fit_vals, color = 'Black')
# plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
# plt.xlabel("Simulation value of Rho")
# plt.ylabel("Error for Rho Estimation")
# # plt.ylim(0.000000001, 0.1)
# plt.xscale('log')
# plt.tight_layout()
# plt.savefig(outDir + "autocorr_comparedEstimates_scatter_5k_error.jpg")
# plt.close()

# estimate_df['corrected'] = estimate_df['Est_Rho'] - estimate_df

#first plot

# sns.set(rc={'figure.figsize':(10,5)}, font_scale = 1)


# ax = sns.stripplot(x = 'Sim_Rho', y = 'C1', data = estimate_df, jitter = True,
#     order = ["0.001", "2e-04", "1e-04", "2e-05", "1e-05", "2e-06"])
# plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
# plt.xlabel("Simulation Value of Rho")
# plt.ylabel("Estimated Value of C1")
# plt.ylim(0,2)
# plt.tight_layout()
# plt.savefig(outDir + "autocorr_C1_5k.jpg")
# plt.close()

# ax = sns.stripplot(x = 'Sim_Rho', y = 'C0', data = estimate_df, jitter = True,
#     order = ["0.001", "2e-04", "1e-04", "2e-05", "1e-05", "2e-06"])
# plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
# plt.xlabel("Simulation Value of Rho")
# plt.ylabel("Estimated Value of C0")
# plt.tight_layout()
# plt.savefig(outDir + "autocorr_C0_5k.jpg")
# plt.close()

# ax = sns.stripplot(x = 'Sim_Rho', y = 'C2', data = estimate_df, jitter = True,
#     order = ["0.001", "2e-04", "1e-04", "2e-05", "1e-05", "2e-06"])
# plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
# plt.xlabel("Simulation Value of Rho")
# plt.ylabel("Estimated Value of C2")
# plt.tight_layout()
# plt.savefig(outDir + "autocorr_C2_5k.jpg")
# plt.close()