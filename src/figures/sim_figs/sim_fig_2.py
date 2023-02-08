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
from sklearn.metrics import mean_squared_error
from matplotlib import rcParams

DIST_TIME_MAX = 50000

NUM_BOOTSTRAPS = 1000
# NUM_BOOTSTRAPS = 10
NUM_REPS = 200
NUM_GROUPS = 100
NUM_SHUFFLES = 500
# NUM_SHUFFLES = 10


#Today I am going to prototype an analysis to determine how well we can
#discriminate between two different recombination rates that are a fixed
#distance apart.

dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_10_03_neutral/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/paper/sim_fig_2/'

dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_10_03_neutral/'
outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/paper/sim_fig_2/'

#First I am going to read the data and randomly pair simulations
all_stat_dfs = []

#loop through each of the dataframes for the separate simulations
for curr_data in os.listdir(dataDir):
    #only get the data directories, not hidden files
    if curr_data[0] == '.':
        continue

    #get the information for the current run
    run_info = curr_data.split('_')
    sim_rho = run_info[1]
    sim_rho = sim_rho[3:]
    rep = run_info[-1]
    if int(rep[3:]) >= NUM_REPS:
        continue

    #get the dataframe for the current run
    d_ratio_file = dataDir + curr_data + "/linkage/d_ratio"
    stat_df = pd.read_pickle(d_ratio_file)
    stat_df['rep'] = int(rep[3:])
    stat_df['Sim_Rho'] = sim_rho
    all_stat_dfs.append(stat_df)
all_stat_dfs = pd.concat(all_stat_dfs)
print(min(all_stat_dfs['d_i']))

#Filter the groups so the first d' value is greater than 0.2
all_stat_dfs = all_stat_dfs[all_stat_dfs['d_i'] > 0.2]

#Randomly divide the reps into 10 groups
rep_groups = np.array(range(0, NUM_REPS))
np.random.shuffle(rep_groups)
rep_groups = np.array_split(rep_groups, NUM_GROUPS)


#Make a dictionary to label each group
group_dict = {}
for i in range(len(rep_groups)):
    for j in rep_groups[i]:
        group_dict[j] = i
        
group_labels = [group_dict[x] for x in all_stat_dfs['rep']]
all_stat_dfs['iter_group'] = group_labels

print("Loaded the data Successfully")

########################## Estimating recombination rates #####################
#loop through each of the distance cutoffs
all_stat_dfs = all_stat_dfs[all_stat_dfs['Dist_X_Time'] <= DIST_TIME_MAX]
all_estimate_df = []
all_fit_df = []

#loop through each rho value
for curr_rho in all_stat_dfs['Sim_Rho'].unique():
    curr_rho_stat = all_stat_dfs[all_stat_dfs['Sim_Rho'] == curr_rho]

    #loop through each iter group
    for curr_iteration in range(0, NUM_GROUPS):
        #get the data for the current rho and iteration
        curr_stat_df = curr_rho_stat[curr_rho_stat['iter_group'] == curr_iteration]

        #Get the current estimate
        lower_fit, upper_fit, estimate_df = plne.bootstrap_rho(curr_stat_df,
                                                             NUM_BOOTSTRAPS)
        estimate_df['Group'] = curr_iteration
        estimate_df['Sim_Rho'] = curr_rho
        all_estimate_df.append(estimate_df)


all_estimate_df = pd.concat(all_estimate_df, ignore_index=True)

all_conf_ints = []

#Calculate the confidence intervals for the estimates with low + high VL
for name, group in all_estimate_df.groupby(['Group', 'Sim_Rho']):
    lower_conf = np.quantile(group['Estimated_Rho'], 0.025)
    mid_est = np.quantile(group['Estimated_Rho'], 0.5)
    upper_conf = np.quantile(group['Estimated_Rho'], 0.975)

    #append the confidence interval
    all_conf_ints.append([lower_conf, mid_est, upper_conf, name[0], name[1]])

all_conf_ints = pd.DataFrame(all_conf_ints, 
    columns=['lower_conf', 'est_rho', 'upper_conf', 'Group', 'Sim_Rho'])

# print(all_conf_ints)


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
for entry in all_conf_ints['Sim_Rho']:
    intRhoList.append(rhoDict[entry])
    newStringRho.append(rho_dict_fix_strings[entry])
all_conf_ints['Sim_float_rho'] = intRhoList
all_conf_ints['Sim_Rho'] = newStringRho

########################## Plotting the accuracy of the estimates #############
fig = sns.stripplot(x = 'Sim_Rho', y = 'est_rho', data = all_conf_ints, 
    jitter = True, color = 'k', s = 1, alpha = 0.3,
    order = [r"$2\times10^{-6}$", r"$10^{-5}$", r"$2\times10^{-5}$", r"$10^{-4}$", r"$2\times10^{-4}$", r"$10^{-3}$"])

# distance across the "X" or "Y" stipplot column to span, in this case 40%
label_width = 0.4
plt.savefig(outDir + "Accuracy.png")
            
for tick, text in zip(fig.get_xticks(), fig.get_xticklabels()):
    sample_name = text.get_text()  # "X" or "Y"

    #get the float value of rho corresponding with the tick
    rho_val = all_conf_ints[all_conf_ints['Sim_Rho'] == sample_name]
    rho_val = rho_val['Sim_float_rho'].unique()[0]

    # plot horizontal lines across the column, centered on the tick
    fig.plot([tick-label_width/2, tick+label_width/2], [rho_val, rho_val],
            lw=2, color='k')
plt.savefig(outDir + "Accuracy.png")
plt.close()


fig = sns.stripplot(x = 'Sim_Rho', y = 'est_rho', data = all_conf_ints, 
    jitter = True, color = 'k', s = 1, alpha = 0.3,
    order = [r"$2\times10^{-6}$", r"$10^{-5}$", r"$2\times10^{-5}$", r"$10^{-4}$", r"$2\times10^{-4}$", r"$10^{-3}$"])

# distance across the "X" or "Y" stipplot column to span, in this case 40%
label_width = 0.4
plt.savefig(outDir + "Accuracy_logy.png")
            
for tick, text in zip(fig.get_xticks(), fig.get_xticklabels()):
    sample_name = text.get_text()  # "X" or "Y"

    #get the float value of rho corresponding with the tick
    rho_val = all_conf_ints[all_conf_ints['Sim_Rho'] == sample_name]
    rho_val = rho_val['Sim_float_rho'].unique()[0]

    # plot horizontal lines across the column, centered on the tick
    fig.plot([tick-label_width/2, tick+label_width/2], [rho_val, rho_val],
            lw=2, color='k')
plt.yscale('log')
plt.savefig(outDir + "Accuracy_logy.png")
plt.close()
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
    curr_pair_diff = group['pair_diff'].unique()[0]
    curr_str_rho_1 = group['str_rho_1'].unique()[0]
    curr_str_rho_2 = group['str_rho_2'].unique()[0]


    #calculate the proportion correct
    correct_order = group[group['order_correct'] == True]
    prop_correct = len(correct_order)/len(group)
    correct_no_overlap = correct_order[correct_order['overlap'] == False]
    prop_correct_no_overlap = len(correct_no_overlap)/len(group)

    disc_results.append([name[0], curr_str_rho_1, name[1], curr_str_rho_2,
             prop_correct, prop_correct_no_overlap, curr_pair_diff])

disc_results = pd.DataFrame(disc_results, columns=['rho_1', 'str_rho_1',
         'rho_2', 'str_rho_2', 'prop_correct', 'prop_correct_no_overlap',
          'pair_diff'])

disc_results['pair_diff'] = disc_results['pair_diff'].apply(lambda x: round(x))
disc_results['pair_diff'] = disc_results['pair_diff'].astype('category')
########################## Plotting our Results ###############################          
print(disc_results['pair_diff'].unique())
rcParams.update({'font.size': 8, 'figure.figsize':(6.5, 2.5)})
sns.set_style("white")

#coolwarm
fig, ax = plt.subplots(1, 2)
sns.stripplot(x = 'str_rho_1', y = 'prop_correct', data = disc_results, 
    jitter = 0.2, s = 3, hue = 'pair_diff', ax = ax[0],
    palette=sns.color_palette("coolwarm", n_colors = len(disc_results['pair_diff'].unique())),
    order = [r"$2\times10^{-6}$", r"$10^{-5}$", r"$2\times10^{-5}$", r"$10^{-4}$", r"$2\times10^{-4}$", r"$10^{-3}$"])
ax[0].axhline(0.5, linestyle= "dashed", color = "black")
ax[0].set_ylim(0,1.1)
ax[0].set_ylabel('Proportion Ordered Correctly')
ax[0].set_xlabel(r'Simulated Recombination Rate ($\rho$)')
ax[0].get_legend().remove()
    

sns.stripplot(x = 'str_rho_1', y = 'prop_correct_no_overlap', data = disc_results, 
    jitter = 0.2, s = 5, hue = 'pair_diff', ax = ax[1],
    palette=sns.color_palette("coolwarm", n_colors = len(disc_results['pair_diff'].unique())),
    order = [r"$2\times10^{-6}$", r"$10^{-5}$", r"$2\times10^{-5}$", r"$10^{-4}$", r"$2\times10^{-4}$", r"$10^{-3}$"])
ax[1].set_ylim(0,1.1)
ax[1].set_ylabel('Proportion with Non-Overlapping \n Confidence Intervals')
ax[1].set_xlabel(r'Simulated Recombination Rate ($\rho$)')
ax[1].legend(title = r'$\rho$ Ratio')
ax[1].axhline(0.05, linestyle= "dashed", color = "black")

for curr_dot in ax[1].legend_.legendHandles:
    curr_dot._sizes = [12]
plt.tight_layout()
plt.savefig(outDir + 'discrimination_conf_refit.png', dpi = 300)
plt.close()