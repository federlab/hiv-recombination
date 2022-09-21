import sys
from turtle import filling
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
from sklearn.metrics import mean_squared_error
from matplotlib import rcParams

THRESHOLD = 0.2
SEG_LOCI = [50, 100, 250, 498]

dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_02_24/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/paper/fig2/'

estimate_df = [] 

#loop through each of the dataframes for the separeate simulations
for curr_data in os.listdir(dataDir):
    #only get the data directories, not hidden files
    if curr_data[0] == '.':
        continue

    #get the information for the current run
    run_info = curr_data.split('_')
    sim_rho = run_info[1]
    sim_rho = sim_rho[3:]
    rep = run_info[-1]

    #get the dataframe for the current run
    d_ratio_file = dataDir + curr_data + "/linkage/d_ratio"
    stat_df = pd.read_pickle(d_ratio_file)

    #loop through the different numbers of segregating loci
    for curr_sample_size in SEG_LOCI:
        #make a list of all the segregating loci
        all_seg_loc_1 = set(stat_df["Locus_1"].unique())
        all_seg_loc_2 = set(stat_df["Locus_2"].unique())
        all_seg_loc = all_seg_loc_1.union(all_seg_loc_2)
        all_seg_loc = np.array(list(all_seg_loc))

        #sample a given number of segregating loci
        seg_loc_sample = np.random.choice(all_seg_loc, size = curr_sample_size,
                             replace = False)
        seg_loc_sample = set(seg_loc_sample)

        #get only the autocorrelation of the chosen loci
        stat_df_sample = stat_df[stat_df["Locus_1"].isin(seg_loc_sample)]
        stat_df_sample = \
            stat_df_sample[stat_df_sample["Locus_2"].isin(seg_loc_sample)]

        #get the estimate and fit for the current dataset and sample size
        x_vals = stat_df_sample['Dist_X_Time'].unique()
        coeffs, fit_dat = optimize.curve_fit(plne.neher_leitner, 
            stat_df_sample['Dist_X_Time'], stat_df_sample['d_ratio'],
            p0 = [0, 0.26, .0000439], maxfev = 10000)
        fit_vals = [plne.neher_leitner(x, coeffs[0], coeffs[1], coeffs[2])
                    for x in x_vals]

        #Bin the d' ratios so they are easier to view on the plots
        binned_rat, binedges, bin_nums = binned_statistic(
            stat_df_sample['Dist_X_Time'].to_numpy(), 
            stat_df_sample['d_ratio'].to_numpy(), bins = 100)

        estimate_df.append([coeffs[0], coeffs[1], coeffs[2], 
            coeffs[1] * coeffs[2], curr_data, sim_rho, curr_sample_size, 
            curr_data])  

estimate_df = pd.DataFrame(estimate_df, columns=["C0", "C1", "C2",
                     "Est_Rho", 'Dataset', 'Sim_Rho', 'Sample Size' , 'data'] )

############################# Plotting Estimate Accuracy ######################
# #Plot our estimates against each other 

#make the rho values ints rather than strings
rhoDict = {"0.001" : 0.001,
            "1e-04" : 0.0001,
            "2e-04" : 0.0002,
            "1e-05" : 0.00001,
            "2e-05" : 0.00002,
            "2e-06" : 0.000002}
rho_dict_fix_strings = { "0.001" : r"$10^{-3}$",
                        "1e-04" : r"$10^{-4}$",
                        "2e-04" : r"$2\times10^{-4}$",
                        "1e-05" : r"$10^{-5}$",
                        "2e-05" : r"$2\times10^{-5}$",
                        "2e-06" : r"$2\times10^{-6}$"}

#redo the labeling on the rho values from what was used in the simulation names
intRhoList = []
newStringRho = []
for entry in estimate_df['Sim_Rho']:
    intRhoList.append(rhoDict[entry])
    newStringRho.append(rho_dict_fix_strings[entry])
estimate_df['Sim_int_rho'] = intRhoList
estimate_df['Sim_Rho'] = newStringRho
print(np.unique(estimate_df['Sim_Rho']))
print(estimate_df[estimate_df['Sim_int_rho'] == 0.000002])

rcParams['mathtext.fontset'] = 'custom'
rcParams['mathtext.rm'] = 'DejaVu Sans'
rcParams['mathtext.it'] = 'DejaVu Sans:italic'


#plot the estimates to show how accurate they are
plot_subset = estimate_df[estimate_df['Sample Size'] == 100]
sns.set(rc={'figure.figsize':(20,10)}, font_scale = 2, font = '')
fig, axes = plt.subplots(1, 2)
sns.stripplot(x = 'Sim_Rho', y = 'Est_Rho', data = plot_subset, 
    jitter = True, color = 'k', s = 8, ax = axes[0],
    order = [r"$2\times10^{-6}$", r"$10^{-5}$", r"$2\times10^{-5}$", r"$10^{-4}$", r"$2\times10^{-4}$", r"$10^{-3}$"])

# distance across the "X" or "Y" stipplot column to span, in this case 40%
label_width = 0.4

            
for tick, text in zip(axes[0].get_xticks(), axes[0].get_xticklabels()):
    sample_name = text.get_text()  # "X" or "Y"

    #get the float value of rho corresponding with the tick
    rho_val = estimate_df[estimate_df['Sim_Rho'] == sample_name]
    rho_val = rho_val['Sim_int_rho'].unique()[0]

    # plot horizontal lines across the column, centered on the tick
    axes[0].plot([tick-label_width/2, tick+label_width/2], [rho_val, rho_val],
            lw=2, color='k')


axes[0].set_xlabel(r'Simulation Value of $\rho$')
axes[0].set_ylabel(r'Estimated Value of $\rho$')
axes[0].set_ylim(0.000001, 0.01)
axes[0].set_yscale('log')


##################### Plotting the MSE for each sample size ###################
#First we need to group by the sample size and rho
#Then we can calculate the mean squared error
grouped_ests = estimate_df.groupby(['Sample Size', 'Sim_int_rho'])
group_MSE = []
for name,group in grouped_ests:
    truth = [name[1] for x in range(len(group))]
    mse = np.sqrt(mean_squared_error(group['Est_Rho'], truth))/np.mean(group['Est_Rho'])
    string_rho = group['Sim_Rho'].unique()[0]
    group_MSE.append([name[0], name[1], mse, string_rho])

group_MSE = pd.DataFrame(group_MSE, 
    columns=['Sample Size', 'Sim_int_rho', 'NRMSE', 'Sim_Rho'])

sns.lineplot(x = 'Sim_int_rho', y = 'NRMSE', 
    data = group_MSE, hue = 'Sample Size', ax = axes[1],
    palette=sns.color_palette("colorblind", n_colors=4))
axes[1].set_xscale('log')
axes[1].set_ylabel('Normalized RMSE')
axes[1].set_xlabel(r'Simulation Value of $\rho$')
plt.legend(title = 'Segregating Loci')
plt.tight_layout()

plt.savefig(outDir + "figure_2.jpg")
plt.close()

    