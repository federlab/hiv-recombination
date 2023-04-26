import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic
from matplotlib import rcParams
import plot_neher as plne
from scipy.stats import linregress

THRESHOLD = 0.2
DIST_TIME_MAX = 50000
NUM_REPS = 200
NUM_GROUPS = 40
NUM_BOOTSTRAPS = 1000

#Today I am going to calculate R^2 for each of the fits and maybe fit a line when R^2 is bad

dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/paper/simulated_estimates/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/paper/supp_line/'

# dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/paper/simulated_estimates/'
# outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/paper/supp_line/'

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
all_stat_dfs = pd.read_pickle(dataDir + "all_stat_dfs_"+ \
                                        str(NUM_BOOTSTRAPS) + "_" + \
                                            str(NUM_GROUPS) + ".pkl")

all_stat_dfs = label_rho(all_stat_dfs)


########################## Estimating recombination rates #####################
#loop through each of the distance cutoffs
all_stat_dfs = all_stat_dfs[all_stat_dfs['Dist_X_Time'] <= DIST_TIME_MAX]
all_estimate_df = []
all_best_fits = []

#Choose a group to plot and initialize the plot
example_group = np.random.choice(range(0, NUM_GROUPS))
#plot the estimates to show how accurate they are
rcParams.update({'figure.figsize':(6.5,8), 'font.size': 8})
fig, axs = plt.subplots(3, 2, sharex = True, sharey = True)
iter_ax = axs.flatten()

#loop through each rho value
for i in range(len(all_stat_dfs['Sim_Rho'].unique())):
    curr_rho = [r"$2\times10^{-6}$", r"$10^{-5}$", r"$2\times10^{-5}$", r"$10^{-4}$", r"$2\times10^{-4}$" , r"$10^{-3}$"][i]
    curr_ax = iter_ax[i]
    curr_rho_stat = all_stat_dfs[all_stat_dfs['Sim_Rho'] == curr_rho]
    curr_float_rho = curr_rho_stat['Sim_float_rho'].unique()[0]

    #loop through each iter group
    for curr_iteration in range(0, NUM_GROUPS):
        #get the data for the current rho and iteration
        curr_stat_df = curr_rho_stat[curr_rho_stat['iter_group'] == curr_iteration]

        #Now downsample the group to the coverage of the lowest bin
        #Bin the d' ratios so they are easier to view on the plots
        bin_count, binedges, bin_nums = binned_statistic(
        curr_stat_df['Dist_X_Time'].to_numpy(), 
        curr_stat_df['d_ratio'].to_numpy(), bins = 100, statistic= 'count')


        #get the estimate and fit for the current dataset and sample size
        x_vals = curr_stat_df['Dist_X_Time'].unique()

        
        line_fit = linregress(curr_stat_df['Dist_X_Time'], curr_stat_df['d_ratio'])
        mid_dist = curr_stat_df[curr_stat_df['Dist_X_Time'] <= 10000]
        mid_dist_fit = linregress(mid_dist['Dist_X_Time'], mid_dist['d_ratio'])
        short_dist = curr_stat_df[curr_stat_df['Dist_X_Time'] <= 2000]
        short_dist_fit = linregress(short_dist['Dist_X_Time'], short_dist['d_ratio'])
        super_short_dist = curr_stat_df[curr_stat_df['Dist_X_Time'] <= 1000]
        super_short_dist_fit = linregress(super_short_dist['Dist_X_Time'], super_short_dist['d_ratio'])
        shortest_dist = curr_stat_df[curr_stat_df['Dist_X_Time'] <= 500]
        shortest_dist_fit = linregress(shortest_dist['Dist_X_Time'], shortest_dist['d_ratio'])

        #put all of the fits into a dataframe
        best_fit_df = []
        best_fit_df.append([shortest_dist_fit.slope, shortest_dist_fit.rvalue**2, 500, shortest_dist_fit.intercept])
        best_fit_df.append([super_short_dist_fit.slope, short_dist_fit.rvalue**2, 1000, super_short_dist_fit.intercept])
        best_fit_df.append([short_dist_fit.slope, short_dist_fit.rvalue**2, 2000, short_dist_fit.intercept])
        best_fit_df.append([line_fit.slope, line_fit.rvalue**2, 50000, line_fit.intercept])
        best_fit_df.append([mid_dist_fit.slope, mid_dist_fit.rvalue**2, 10000, mid_dist_fit.intercept])

        best_fit_df = pd.DataFrame(best_fit_df, columns = ['Estimated_Rho', 'Fit_R2', 'Fit Length', 'Intercept'])
        best_fit_df['Sim_Rho'] = curr_rho
        all_best_fits.append(best_fit_df)

        #Plot the fit for an example group
        if curr_iteration == example_group:
            plot_df = []
            for ind, row in best_fit_df.iterrows():
                fit_vals = row['Estimated_Rho'] * x_vals + row['Intercept']
                fit_vals = pd.DataFrame(list(zip(fit_vals, x_vals)), columns = ['Fit', 'Dist_X_Time'])
                fit_vals['Fit Length'] = row['Fit Length']
                plot_df.append(fit_vals)
            plot_df = pd.concat(plot_df, ignore_index = True)

            #Bin the d' ratios so they are easier to view on the plots
            binned_rat, binedges, bin_nums = binned_statistic(
                curr_stat_df['Dist_X_Time'].to_numpy(), 
                curr_stat_df['d_ratio'].to_numpy(), bins = 100)

            sns.lineplot(y = binned_rat, x = binedges[1:], color = 'black', ax = curr_ax)
            sns.lineplot(x = 'Dist_X_Time', y = 'Fit', hue = 'Fit Length', data = plot_df, ax = curr_ax,
                        palette = sns.color_palette("coolwarm", n_colors = best_fit_df.shape[0]))
            curr_ax.set_title(r'$\rho$ = ' + str(curr_rho), fontsize = 8)
            curr_ax.set_ylim(-0.1,1.5)
            curr_ax.set_ylabel('D\' Ratio')
            curr_ax.set_xlabel(r'Distance $\cdot$ Time (bp $\cdot$ generations)')

            #Remove the legends except for the last one
            if i != len(all_stat_dfs['Sim_Rho'].unique()) - 1:
                curr_ax.get_legend().remove()

plt.savefig(outDir + 'resampled_fits.png', dpi = 300)
plt.close()


all_best_fits = pd.concat(all_best_fits, ignore_index=True)


#To get the comparison, we need to convert the rho values to floats
rhoDict = {r"$10^{-3}$" : 0.001,
            r"$10^{-4}$" : 0.0001,
            r"$2\times10^{-4}$" : 0.0002,
            r"$10^{-5}$" : 0.00001,
            r"$2\times10^{-5}$" : 0.00002,
            r"$2\times10^{-6}$" : 0.000002}


#redo the labeling on the rho values from what was used in the simulation names
intRhoList = []
newStringRho = []
for entry in all_best_fits['Sim_Rho']:
    intRhoList.append(rhoDict[entry])
all_best_fits['Sim_float_rho'] = intRhoList
print(len(all_best_fits))

negative_ests = all_best_fits[all_best_fits['Estimated_Rho'] < 0]
print(negative_ests['Fit Length'].value_counts())
# print(len(negative_ests))
# print(negative_ests['Sim_float_rho'].unique())
# print(negative_ests['Fit Length'].unique())

########################## Plotting the accuracy of the estimates #############
#plot the estimates to show how accurate they are
rcParams.update({'figure.figsize':(3.5,3), 'font.size': 8})

#make a 

fig = sns.stripplot(x = 'Sim_Rho', y = 'Estimated_Rho', data = all_best_fits,
    jitter = True, hue = 'Fit Length', s = 3, alpha = 1,
    order = [r"$2\times10^{-6}$", r"$10^{-5}$", r"$2\times10^{-5}$", r"$10^{-4}$", r"$2\times10^{-4}$" , r"$10^{-3}$"],
    palette=sns.color_palette("coolwarm", n_colors=len(all_best_fits['Fit Length'].unique())))

# distance across the "X" or "Y" stipplot column to span, in this case 40%
label_width = 0.4

#need to save the fig to force the axis to be drawn
plt.savefig(outDir + "line_fits_both.png", dpi = 300)
            
for tick, text in zip(fig.get_xticks(), fig.get_xticklabels()):
    sample_name = text.get_text()  # "X" or "Y"

    #get the float value of rho corresponding with the tick
    rho_val = all_best_fits[all_best_fits['Sim_Rho'] == sample_name]
    rho_val = rho_val['Sim_float_rho'].unique()[0]

    # plot horizontal lines across the column, centered on the tick
    fig.plot([tick-label_width/2, tick+label_width/2], [rho_val, rho_val],
            lw=1, color='k', zorder = 3)

fig.set_ylabel(r'Estimated Recombination Rate ($\hat{\rho}$)')
fig.set_xlabel(r'Simulated Recombination Rate ($\rho$)')
fig.set_yscale('log')
plt.legend(loc = 'lower right')
plt.tight_layout()
plt.savefig(outDir + "resampled_line_fits_both_" + str(NUM_GROUPS) + ".png", dpi = 300)
plt.close()