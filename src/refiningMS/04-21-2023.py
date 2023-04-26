import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import estimation_util as est_util
from matplotlib import rcParams
from matplotlib.lines import Line2D
from scipy.stats import binned_statistic


DIST_TIME_MAX = 50000
RATIOS_PER_GROUP = 25000

NUM_BOOTSTRAPS = 1000
NUM_REPS = 200
NUM_GROUPS = 40
NUM_SHUFFLES = 500


#Today I am workshopping a 4 panel figure of the simulated data
#It will include the D' ratio curves, accuracy analysis, and both 
#discrimination analyses

dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/paper/simulated_estimates/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_10_03_neutral/04-21-2023/'

# dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/paper/simulated_estimates/'
# outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_10_03_neutral/04-21-2023/'
###############################################################################
#################################### Helper Functions #########################
###############################################################################

def label_rho(my_dataframe):
    #This function takes a dataframe and labels the rho values with floats
    #instead of strings so that they can be used as a continuous variable
    #To get the comparison, we need to convert the rho values to floats
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
    for entry in my_dataframe['Sim_Rho']:
        intRhoList.append(rhoDict[entry])
        newStringRho.append(rho_dict_fix_strings[entry])
    my_dataframe['Sim_float_rho'] = intRhoList
    my_dataframe['Sim_Rho'] = newStringRho

    return my_dataframe
###############################################################################
###############################################################################
###############################################################################
#Load the dataframes that will be plotted
all_conf_ints = pd.read_pickle(dataDir + "all_conf_ints_"+ \
                                        str(NUM_BOOTSTRAPS) + "_" + \
                                            str(NUM_GROUPS) + ".pkl")
all_stat_dfs = pd.read_pickle(dataDir + "all_stat_dfs_"+ \
                                        str(NUM_BOOTSTRAPS) + "_" + \
                                            str(NUM_GROUPS) + ".pkl")

all_conf_ints = label_rho(all_conf_ints)
all_stat_dfs = label_rho(all_stat_dfs)

######################## Bin the values for the D' Ratio Curves ###############
all_rho_bins = []

#Bin the points for each rho value
for curr_rho in all_stat_dfs['Sim_Rho'].unique():
    curr_stat_df = all_stat_dfs[all_stat_dfs['Sim_Rho'] == curr_rho]

    #Bin the d' ratios so they are easier to view on the plots
    binned_rat, binedges, bin_nums = binned_statistic(
        curr_stat_df['Dist_X_Time'].to_numpy(), 
        curr_stat_df['d_ratio'].to_numpy(), bins = 100)
    
    curr_rho_bins = pd.DataFrame({'Dist_X_Time': binedges[1:], 'D_Ratio': binned_rat})
    curr_rho_bins['Sim_Rho'] = curr_rho
    curr_rho_bins['Sim_float_rho'] = curr_stat_df['Sim_float_rho'].unique()[0]

    all_rho_bins.append(curr_rho_bins)

all_rho_bins = pd.concat(all_rho_bins)

######################## Running Discrimination Analysis ######################

#We need to loop through all pairs of rates
possible_rhos = all_conf_ints['Sim_float_rho'].unique()
possible_rhos = np.sort(possible_rhos)

pair_results_df = []

#loop through all but the highest rho value
for i in range(len(possible_rhos -1)):
    curr_rho_1 = possible_rhos[i]
    rho_1_df = all_conf_ints[all_conf_ints['Sim_float_rho'] == curr_rho_1]
    curr_str_rho_1 = rho_1_df['Sim_Rho'].unique()[0]

    for j in range(i, len(possible_rhos)):
        curr_rho_2 = possible_rhos[j]
        rho_2_df = all_conf_ints[all_conf_ints['Sim_float_rho'] == curr_rho_2]
        curr_str_rho_2 = rho_2_df['Sim_Rho'].unique()[0]

        pair_diff = curr_rho_2/curr_rho_1
        if pair_diff >= 11:
            continue

        for curr_shuffle in range(NUM_SHUFFLES):
            #randomly pair the simulations
            rep_list_1 = rho_1_df['Group'].unique()
            np.random.shuffle(rep_list_1)
            rep_list_2 = rho_2_df['Group'].unique()
            np.random.shuffle(rep_list_2)
            assert len(rep_list_1) == len(rep_list_2), "The number of replicates for each rho value is not the same"
            rep_pairs = [(rep_list_1[i], rep_list_2[i]) for i in range(len(rep_list_1))]
            
            #loop through each pair and check if the confidence intervals overlap
            for curr_pair in rep_pairs:
                rho_1_group = rho_1_df[rho_1_df['Group'] == curr_pair[0]]
                rho_2_group = rho_2_df[rho_2_df['Group'] == curr_pair[1]]

                #check if the estimates are in the correct order
                if rho_1_group['est_rho'].values[0] > rho_2_group['est_rho'].values[0]:
                    order_correct = False
                else:
                    order_correct = True  
                
                #check if the confidence intervals overlap
                lower_conf_1 = rho_1_group['lower_conf'].values[0]
                upper_conf_1 = rho_1_group['upper_conf'].values[0]
                lower_conf_2 = rho_2_group['lower_conf'].values[0]
                upper_conf_2 = rho_2_group['upper_conf'].values[0]
                
                if (lower_conf_1 < lower_conf_2) and (lower_conf_2 < upper_conf_1):
                    overlap = True
                elif (lower_conf_1 < upper_conf_2) and (upper_conf_2 < upper_conf_1):
                    overlap = True
                else:
                    overlap = False

                #Don't compare reps to themselves
                if curr_pair[0] == curr_pair[1] and round(pair_diff) ==1:
                    continue

                #append the results
                pair_results_df.append([curr_rho_1, curr_str_rho_1, curr_rho_2, 
                            curr_str_rho_2, order_correct, overlap, pair_diff])

pair_results_df = pd.DataFrame(pair_results_df, columns=['rho_1', 'str_rho_1',
        'rho_2', 'str_rho_2', 'order_correct', 'overlap', 'pair_diff'])

#Now calculate the proportion correct 
disc_results = []
pair_results_df = pair_results_df.groupby(['rho_1', 'rho_2'])

for name, group in pair_results_df:
    ordered_rhos = [name[0], name[1]]
    ordered_rhos.sort()
    curr_pair_diff = group['pair_diff'].unique()[0]

    #Get the correct string labels for the paired rho values
    reverse_rhoDict = {0.001 : r"$10^{-3}$",
                0.0001 : r"$10^{-4}$",
                0.0002 : r"$2\times10^{-4}$",
                0.00001 : r"$10^{-5}$",
                0.00002 : r"$2\times10^{-5}$",
                0.000002 : r"$2\times10^{-6}$"}
    curr_str_rho_1 = reverse_rhoDict[ordered_rhos[0]]
    curr_str_rho_2 = reverse_rhoDict[ordered_rhos[1]]

    #calculate the proportion correct
    correct_order = group[group['order_correct'] == True]
    prop_correct = len(correct_order)/len(group)
    percent_correct = round(prop_correct * 100, 0)
    correct_no_overlap = correct_order[correct_order['overlap'] == False]
    prop_correct_no_overlap = len(correct_no_overlap)/len(group)
    percent_correct_no_overlap = round(prop_correct_no_overlap * 100,0)

    disc_results.append([ordered_rhos[0], curr_str_rho_1, ordered_rhos[1], curr_str_rho_2,
             prop_correct, percent_correct, prop_correct_no_overlap, 
             percent_correct_no_overlap, curr_pair_diff])

disc_results = pd.DataFrame(disc_results, columns=['rho_1', 'str_rho_1',
         'rho_2', 'str_rho_2', 'prop_correct', 'percent_correct',
         'prop_correct_no_overlap', 'percent_correct_no_overlap',
          'pair_diff'])


disc_results['pair_diff'] = disc_results['pair_diff'].apply(lambda x: round(x))
disc_results['pair_diff'] = disc_results['pair_diff'].astype('category')

label_dict = {
    1 : r"$\rho_1$",
    2 : r"$2 \times \rho_1$",
    5 : r"$5\times \rho_1$",
    10 : r"$10 \times \rho_1$",
}

label_dict = {
    1 : r"$\rho$ vs $\rho$",
    2 : r"$\rho$ vs $2\rho$",
    5 : r"$\rho$ vs $5\rho$",
    10 : r"$\rho$ vs $10\rho$",
}

disc_results['pair_diff'] = disc_results['pair_diff'].map(label_dict)

######################### Configure Plot Settings #############################
#plot the estimates to show how accurate they are
params = {'figure.figsize':(6, 4.5), 'axes.labelsize': 6,'axes.titlesize':6,  
          'legend.fontsize': 6, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
          'legend.title_fontsize': 6}
rcParams.update(params)

#Make a 4 panel Figure
fig, axs = plt.subplots(2, 2)
plt.subplots_adjust(hspace = 0.4, wspace = 0.35)

######################### Panel A  D' Ratio ###################################
sns.lineplot(data=all_rho_bins, x='Dist_X_Time', y='D_Ratio', 
             hue='Sim_float_rho', linewidth = 1, ax = axs[0,0],
             palette=sns.color_palette("viridis", 
             n_colors=len(all_rho_bins['Sim_float_rho'].unique())))

axs[0,0].get_legend().remove()
axs[0,0].set_xlabel(r'Distance $\cdot$ Time (bp $\cdot$ generations)')
axs[0,0].set_ylabel("D\' Ratio")
axs[0,0].set_xlim(0,50000)
# axs[0,0].ticklabel_format(axis="x", style="sci", scilimits=(0,0))

######################### Panel B  Accuracy ###################################
sns.stripplot(x = 'Sim_Rho', y = 'est_rho', hue = 'Sim_float_rho',
                data = all_conf_ints, jitter = True, s = 4,
                palette=sns.color_palette("viridis", 
                n_colors=len(all_rho_bins['Sim_float_rho'].unique())),
                alpha = 0.3, ax = axs[0, 1], order = [r"$2\times10^{-6}$", 
                r"$10^{-5}$", r"$2\times10^{-5}$", r"$10^{-4}$", 
                r"$2\times10^{-4}$", r"$10^{-3}$"])

# distance across the "X" or "Y" stipplot column to span, in this case 40%
label_width = 0.5

# we need to save the figure once here so that we can access the axis labels
plt.savefig(outDir + "figure3_prototype.jpg", dpi = 300)
            
for tick, text in zip(axs[0,1].get_xticks(), axs[0,1].get_xticklabels()):
    sample_name = text.get_text()  # "X" or "Y"

    #get the float value of rho corresponding with the tick
    rho_val = all_conf_ints[all_conf_ints['Sim_Rho'] == sample_name]
    rho_val = rho_val['Sim_float_rho'].unique()[0]
    estimate = np.mean(all_conf_ints[all_conf_ints['Sim_Rho'] == sample_name]['est_rho'])

    # plot horizontal lines across the column, centered on the tick
    axs[0,1].plot([tick-label_width/2, tick+label_width/2], [rho_val, rho_val],
            lw=1, color='r', linestyle = '--')
    axs[0,1].plot([tick-label_width/2, tick+label_width/2], [estimate, estimate],
            lw=1, color='k')


# make a legend showing which line is the mean and which line is the truth
redline = Line2D([0], [0], color='r', linestyle = '--', label = r'True $\rho$')
blackline = Line2D([0], [0], color='k', label = r'Mean Estimate')
axs[0,1].legend(handles = [redline, blackline], loc = 'upper left', frameon = True)


axs[0,1].set_ylabel("Estimated Recombination "+ r'Rate ($\hat{\rho}$)')
axs[0,1].set_xlabel(r'Simulated Recombination Rate ($\rho$)')
axs[0,1].set_yscale('log')

########################## Panel C Ordering ###################################
sns.stripplot(x = 'pair_diff', y = 'percent_correct', data = disc_results, 
    jitter = 0.3, s = 4, hue = 'str_rho_1', ax = axs[1,0],
    palette=sns.color_palette("viridis", n_colors = len(disc_results['str_rho_1'].unique())),)
axs[1,0].set_ylim(-5,105)
axs[1,0].set_ylabel('Percent of Pairs Correctly Ranked')
axs[1,0].set_xlabel("Compared Recombination Rates")
axs[1,0].axhline(50, linestyle= "dashed", color = "black", linewidth = 1)
axs[1,0].get_legend().remove()

# blackline = Line2D([0], [0], color='k', linestyle='--', label = r'Expected Accuracy for $\rho$ vs $\rho$')
# axs[1,0].legend(handles = [blackline], loc = 'lower right', frameon = True)

########################## Panel D Significance ################################
sns.stripplot(x = 'pair_diff', y = 'percent_correct_no_overlap', data = disc_results, 
    jitter = 0.3, s = 4, hue = 'str_rho_1', ax = axs[1,1],
    palette=sns.color_palette("viridis", n_colors = len(disc_results['str_rho_1'].unique())))
axs[1,1].set_ylim(-5,105)
axs[1,1].set_ylabel('Percent of Differences Deemed Significant')
axs[1,1].set_xlabel("Compared Recombination Rates")
axs[1,1].axhline(5, linestyle= "dashed", color = "black", linewidth = 1)
axs[1,1].get_legend().remove()
# blackline = Line2D([0], [0], color='k', linestyle='--', label = r'5% Significance Level')
# axs[1,1].legend(handles = [blackline], loc = 'lower right', frameon = True)


######################### Format the Figure Legend ############################
#Relabel the legend handles so the values are formatted nicely
new_labels = [r"$2\times10^{-6}$", r"$10^{-5}$", r"$2\times10^{-5}$",
               r"$10^{-4}$", r"$2\times10^{-4}$", r"$10^{-3}$" ]

fig.legend(loc='upper center', bbox_to_anchor=(0.5, .97),title = r'Simulated Recombination Rate ($\rho$)', ncol = 6,
           labels = new_labels, handles = axs[0,0].get_legend_handles_labels()[0], frameon = True)

axs[0,0].set_xticks(axs[0,0].get_xticks(), axs[0,0].get_xticklabels(), rotation=30, ha='right')
axs[0,1].set_xticks(axs[0,1].get_xticks(), axs[0,1].get_xticklabels(), rotation=30, ha='right')
axs[1,0].set_xticks(axs[1,0].get_xticks(), axs[1,0].get_xticklabels(), rotation=30, ha='right')
axs[1,1].set_xticks(axs[1,1].get_xticks(), axs[1,1].get_xticklabels(), rotation=30, ha='right')


plt.savefig(outDir + "figure3_prototype.jpg", dpi=300)