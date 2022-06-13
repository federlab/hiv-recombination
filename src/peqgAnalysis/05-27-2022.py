import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import zaniniUtil as zu
from scipy.stats import binned_statistic
import plot_neher as plne
import autocorrelation as autocorr
from scipy import optimize

#################################################### Figure 2 #######################################################
#Today I am quickly plotting the viral loads to add to my poster

dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/'
outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/zanini/poster/'

dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/zanini/poster_peqg/'

available_files = os.listdir(dataDir)

#data frame of viral loads for each patient
viralLoadData = zu.make_viral_load_df(dataDir)

#now plot the viral loads
sns.set(rc={'figure.figsize':(13.5,8.5)}, font_scale = 2)
myplot = sns.lineplot(x = 'Days from infection', y = 'Viral load [virions/ml]', hue = 'Participant', data = viralLoadData,
                     ci = None)
myplot.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
myplot.set(yscale = 'log')
plt.xlabel("Days from Infection")
plt.ylabel("Viral Load [virions/ml]")
plt.xlim(-10, 3000)
plt.tight_layout()
plt.savefig(outDir + "viralLoads", dpi = 300)
plt.close()

quit()
#################################################### Figure 3 #######################################################
THRESHOLD = 0.2
SEG_LOCI = 100

dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_02_24/'

estimate_df = [] 

for curr_data in os.listdir(dataDir):
    #only get the data directories, not hidden files
    if curr_data[0] == '.':
        continue
    run_info = curr_data.split('_')
    sim_rho = run_info[1]
    sim_rho = sim_rho[3:]
    rep = run_info[-1]

    linkage_file = dataDir + curr_data + "/linkage/r2_and_D"
    d_ratio_file = dataDir + curr_data + "/linkage/d_ratio"

    if os.path.exists(d_ratio_file):
        stat_df = pd.read_pickle(d_ratio_file)
    else:
        stat_df = autocorr.calculate_d_ratios(linkage_file, THRESHOLD)
        stat_df.to_pickle(d_ratio_file)
    
    #make a list of all the segregating loci
    all_seg_loc_1 = set(stat_df["Locus_1"].unique())
    all_seg_loc_2 = set(stat_df["Locus_2"].unique())
    all_seg_loc = all_seg_loc_1.union(all_seg_loc_2)
    all_seg_loc = np.array(list(all_seg_loc))

    #sample a given number of segregating loci
    seg_loc_sample = np.random.choice(all_seg_loc, size = SEG_LOCI, replace = False)
    seg_loc_sample = set(seg_loc_sample)

    #get only the autocorrelation of the chosen loci
    stat_df_sample = stat_df[stat_df["Locus_1"].isin(seg_loc_sample)]
    stat_df_sample = stat_df_sample[stat_df_sample["Locus_2"].isin(seg_loc_sample)]

    x_vals = stat_df_sample['Dist_X_Time'].unique()
    coeffs, fit_dat = optimize.curve_fit(plne.neher_leitner, stat_df_sample['Dist_X_Time'], stat_df_sample['d_ratio'], p0 = [0, 0.26, .0000439])
    fit_vals = [plne.neher_leitner(x, coeffs[0], coeffs[1], coeffs[2]) for x in x_vals]
    # sns.scatterplot(x = 'Dist_X_Time', y = 'd_ratio', data = stat_df_sample)
    # sns.lineplot(x = 'Dist_X_Time', y = 'd_ratio', data = stat_df_sample, estimator = np.mean)
    # sns.lineplot(x = 'Dist_X_Time', y = 'd_ratio', data = stat_df_sample, estimator = np.median)
    # plt.savefig(currOut + "/auto_plot.jpg")
    # plt.close()

    estimate_df.append([coeffs[0], coeffs[1], coeffs[2], coeffs[1] * coeffs[2], curr_data, sim_rho, 'curve', curr_data])  

estimate_df = pd.DataFrame(estimate_df, columns=["C0", "C1", "C2", "Est_Rho", 'Dataset', 'Sim_Rho', 'style', 'data'] )
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

sns.set(rc={'figure.figsize':(15,10)}, font_scale = 2)
ax = sns.stripplot(x = 'Sim_Rho', y = 'Est_Rho', data = estimate_df, jitter = True, color = 'k', s = 8,
    order = ["0.001", "2e-04", "1e-04", "2e-05", "1e-05", "2e-06"])

# distance across the "X" or "Y" stipplot column to span, in this case 40%
label_width = 0.4

for tick, text in zip(ax.get_xticks(), ax.get_xticklabels()):
    sample_name = text.get_text()  # "X" or "Y"

    # calculate the median value for all replicates of either X or Y
    rho_val = rhoDict[sample_name]

    # plot horizontal lines across the column, centered on the tick
    plt.plot([tick-label_width/2, tick+label_width/2], [rho_val, rho_val],
            lw=2, color='k')

plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
plt.xlabel(r'Simulation Value of $\rho$')
plt.ylabel(r'Estimated Value of $\rho$')
plt.ylim(0.000001, 0.01)
plt.tight_layout()
plt.yscale('log')
plt.savefig(outDir + "figure_3.jpg", dpi = 300)