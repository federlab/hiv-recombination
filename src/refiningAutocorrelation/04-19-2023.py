import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import estimation_util as est_util
from matplotlib import rcParams
from matplotlib.lines import Line2D


DIST_TIME_MAX = 50000
RATIOS_PER_GROUP = 25000

NUM_BOOTSTRAPS = 1000
NUM_REPS = 200
NUM_GROUPS = 40
NUM_SHUFFLES = 500


#Today I am reworking the simulation figure to use formatting that is more
#clear to people.

dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/paper/simulated_estimates/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_10_03_neutral/04-19-2023/'

# dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/paper/simulated_estimates/'
# outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_10_03_neutral/04-19-2023/'
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

all_conf_ints = pd.read_pickle(dataDir + "all_conf_ints_"+ \
                                        str(NUM_BOOTSTRAPS) + "_" + \
                                            str(NUM_GROUPS) + ".pkl")

all_conf_ints = label_rho(all_conf_ints)

########################## Comparing pairs of groups ##########################

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
    1 : r"$\rho$"+ "\n vs\n " + r"$\rho$",
    2 : r"$\rho$"+ "\n vs\n " + r"$2\rho$",
    5 : r"$\rho$"+ "\n vs\n " + r"$5\rho$",
    10 : r"$\rho$"+ "\n vs\n " + r"$10\rho$",
}

label_dict = {
    1 : r"$\rho_1$",
    2 : r"$2 \times \rho_1$",
    5 : r"$5\times \rho_1$",
    10 : r"$10 \times \rho_1$",
}

disc_results['pair_diff'] = disc_results['pair_diff'].map(label_dict)


########################## Plotting our Results ###############################          


rcParams.update({'font.size': 8, 'figure.figsize':(7, 3)})
fig, ax = plt.subplots(1, 3)

sns.stripplot(x = 'Sim_Rho', y = 'est_rho', data = all_conf_ints, 
    jitter = True, color = 'k', s = 3, alpha = 0.3, ax = ax[0],
    order = [r"$2\times10^{-6}$", r"$10^{-5}$", r"$2\times10^{-5}$", r"$10^{-4}$", r"$2\times10^{-4}$", r"$10^{-3}$"])

# distance across the "X" or "Y" stipplot column to span, in this case 40%
label_width = 0.5

# we need to save the figure once here so that we can access the axis labels
plt.tight_layout()
plt.savefig(outDir + "intermediate.png", dpi = 300)

params = {'figure.figsize':(3, 3), 'axes.labelsize': 8,'axes.titlesize':8,  'legend.fontsize': 8, 'xtick.labelsize': 8, 'ytick.labelsize': 8}
rcParams.update(params)
            
for tick, text in zip(ax[0].get_xticks(), ax[0].get_xticklabels()):
    sample_name = text.get_text()  # "X" or "Y"

    #get the float value of rho corresponding with the tick
    rho_val = all_conf_ints[all_conf_ints['Sim_Rho'] == sample_name]
    rho_val = rho_val['Sim_float_rho'].unique()[0]
    estimate = np.mean(all_conf_ints[all_conf_ints['Sim_Rho'] == sample_name]['est_rho'])

    # plot horizontal lines across the column, centered on the tick
    ax[0].plot([tick-label_width/2, tick+label_width/2], [rho_val, rho_val],
            lw=1, color='r', linestyle = '--')
    ax[0].plot([tick-label_width/2, tick+label_width/2], [estimate, estimate],
            lw=1, color='k')


# make a legend showing which line is the mean and which line is the truth
redline = Line2D([0], [0], color='r', linestyle = '--', label = r'True $\rho$')
blackline = Line2D([0], [0], color='k', label = r'Mean Estimate')
ax[0].legend(handles = [redline, blackline], loc = 'upper left', frameon = True)


ax[0].set_ylabel(r'Estimated Recombination Rate ($\hat{\rho}$)')
ax[0].set_xlabel(r'Simulated Recombination Rate ($\rho$)')
ax[0].set_yscale('log')

sns.stripplot(x = 'pair_diff', y = 'percent_correct', data = disc_results, 
    jitter = 0.5, s = 3, hue = 'str_rho_1', ax = ax[1],
    palette=sns.color_palette("viridis", n_colors = len(disc_results['str_rho_1'].unique())),)
ax[1].set_ylim(-5,105)
ax[1].set_ylabel('Percent of Estimates Correctly Ordered')
ax[1].set_xlabel('Second Recombination Rate ' + r"($\rho_2$)")
ax[1].axhline(50, linestyle= "dashed", color = "black")
ax[1].legend(loc='lower right', borderaxespad=0, title = "First Recomb. \n   Rate " + r"$\rho_1$")

sns.stripplot(x = 'pair_diff', y = 'percent_correct_no_overlap', data = disc_results, 
    jitter = 0.5, s = 3, hue = 'str_rho_1', ax = ax[2],
    palette=sns.color_palette("viridis", n_colors = len(disc_results['str_rho_1'].unique())))
ax[2].set_ylim(-5,105)
ax[2].set_ylabel("Percent of Differences \n Deemed Significant")
ax[2].set_xlabel('Second Recombination Rate ' + r"($\rho_2$)")
ax[2].axhline(5, linestyle= "dashed", color = "black")
ax[2].get_legend().remove()

for curr_dot in ax[1].legend_.legendHandles:
    curr_dot._sizes = [12]

ax[0].set_xticks(ax[0].get_xticks(), ax[0].get_xticklabels(), rotation=45, ha='right')
ax[1].set_xticks(ax[1].get_xticks(), ax[1].get_xticklabels(), rotation=45, ha='right')
ax[2].set_xticks(ax[2].get_xticks(), ax[2].get_xticklabels(), rotation=45, ha='right')

plt.tight_layout()
plt.savefig(outDir + "option2.png", dpi = 300)
plt.close()