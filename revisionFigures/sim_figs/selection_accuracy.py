import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import slimUtil as slim_util
from matplotlib import rcParams
from matplotlib.lines import Line2D
from scipy.stats import binned_statistic



# This figure will be a three panel of the accuracy for the selection
# simulations, vNE simulations, and population substructure
dataDirSel = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/revision/simulated_estimates_selection/'
dataDirSub = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/revision/simulated_estimates_sub/'
dataDirSubInset = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/revision/simulated_linkage_sub/'
dataDirSelInset = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/revision/simulated_linkage_sel/'
dataDirVNE = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/revision/simulated_estimates_vNE/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/revision/sim_figs/'
outFileName = 'selection_accuracy.jpg'

###############################################################################
##################### Setting the figure parameters ###########################
###############################################################################
#plot the estimates to show how accurate they are
params = {'figure.figsize':(7, 2.5), 'axes.labelsize': 6,'axes.titlesize':6,  
          'legend.fontsize': 6, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
          'legend.title_fontsize': 6}
rcParams.update(params)
fig, axs = plt.subplots(1, 3, gridspec_kw={'bottom' : 0.3},
                 sharey= True)
plt.subplots_adjust(wspace = 0.05)

label_width = 0.5
linewidth = 0.75
markersize = 4
markeredgewidth = 1


###############################################################################
########################## Selection Accuracy Panel ###########################
###############################################################################
SELECTION_RHO_LIST = [r"$2\times10^{-6}$",
                        r"$5\times10^{-6}$",
                        r"$10^{-5}$",
                        r"$2\times10^{-5}$",
                        r"$5\times10^{-5}$",
                        r"$10^{-4}$",
                        r"$2\times10^{-4}$",
                        r"$5\times10^{-4}$",
                        r"$10^{-3}$"]
NUM_BOOTSTRAPS_SEL = 1000
NUM_GROUPS_SEL = 10

all_conf_ints_sel = pd.read_pickle(dataDirSel+ "all_conf_ints_"+ \
                                        str(NUM_BOOTSTRAPS_SEL) + "_" + \
                                            str(NUM_GROUPS_SEL) + ".pkl")

#Relabel the Rho values
all_conf_ints_sel = slim_util.label_rho(all_conf_ints_sel)
all_conf_ints_sel = all_conf_ints_sel[
        all_conf_ints_sel['Sim_Rho'].isin(SELECTION_RHO_LIST)]

all_conf_ints_sel['Estimated_Rho'] = all_conf_ints_sel['est_rho']

#Now plot the accuracy
ax0 = axs[0]
sns.stripplot(x = 'Sim_Rho', y = 'Estimated_Rho', data = all_conf_ints_sel,
               jitter = True, s = 2, color = 'k', alpha = 0.3,
                order = SELECTION_RHO_LIST, ax = ax0)

plt.savefig(outDir + outFileName, dpi = 400)
            
for tick, text in zip(ax0.get_xticks(), ax0.get_xticklabels()):
    sample_name = text.get_text()  # "X" or "Y"

    #get the float value of rho corresponding with the tick
    rho_val = all_conf_ints_sel[all_conf_ints_sel['Sim_Rho'] == sample_name]
    rho_val = rho_val['Sim_float_rho'].unique()[0]
    estimate = np.mean(all_conf_ints_sel[
            all_conf_ints_sel['Sim_Rho'] == sample_name]['Estimated_Rho'])

    # plot horizontal lines across the column, centered on the tick
    ax0.plot([tick-label_width/2, tick+label_width/2], [rho_val, rho_val],
            lw= linewidth, color='r', linestyle = 'dotted')
    ax0.plot([tick-label_width/2, tick+label_width/2], [estimate, estimate],
            lw= linewidth, color='k')

redline = Line2D([0], [0], color='r', linestyle = 'dotted',
                  label = r'True $\rho$', linewidth= linewidth)
blackline = Line2D([0], [0], color='k', label = r'Mean Estimate',
                    linewidth= linewidth)    
ax0.set_title('Selection')
ax0.set_ylim(0.0000001, 0.01)
ax0.set_xticklabels(ax0.get_xticklabels(), rotation = 90)
ax0.set_xlabel(r'True Recombination Rate ($\rho$)')
ax0.set_ylabel(r'Estimated Recombination Rate $(\hat{\rho})$')
ax0.set_yscale('log')
ax0.legend(handles = [redline, blackline], loc = 'upper left', frameon = False,
           columnspacing = 0, handletextpad = 0.1) 

###############################################################################
######################### Substructure Accuracy Panel #########################
###############################################################################
SUB_RHO_LIST = [r"$10^{-5}$", r"$10^{-4}$"]
NUM_BOOTSTRAPS_SUB = 1000
NUM_GROUPS_SUB = 5

all_conf_ints_sub = pd.read_pickle(dataDirSub + "all_conf_ints_"+ \
                                        str(NUM_BOOTSTRAPS_SUB) + "_" + \
                                            str(NUM_GROUPS_SUB) + ".pkl")

#Relabel the Rho values and rename the accuracy column
all_conf_ints_sub = slim_util.label_rho(all_conf_ints_sub)
all_conf_ints_sub['Estimated_Rho'] = all_conf_ints_sub['est_rho']

#Make a subset of the selection data to include on this plot
subset_sel = all_conf_ints_sel[all_conf_ints_sel['Sim_Rho'].isin(SUB_RHO_LIST)]
subset_sel = subset_sel[subset_sel['Group'] < NUM_GROUPS_SUB]
subset_sel['m_rate'] = 'Well Mixed'

all_conf_ints_sub = pd.concat(
                        [all_conf_ints_sub, subset_sel], ignore_index=True)
all_conf_ints_sub = slim_util.label_m_rate(all_conf_ints_sub)
all_conf_ints_sub['m_rate'] = pd.Categorical(all_conf_ints_sub['m_rate'],
                                categories = ["Well Mixed", r"$m=10^{-2}$",
                                r"$m=10^{-3}$", r"$m=10^{-4}$", r"$m=10^{-5}$"])

#Now Plot the accuracy
ax1 = axs[1]
sns.stripplot(x = 'Sim_Rho', y = 'Estimated_Rho',
                data = all_conf_ints_sub, jitter = True, s = 2,
                palette = sns.color_palette("YlOrRd", 
                n_colors= len(all_conf_ints_sub['m_rate'].unique())),
                hue = 'm_rate', order = SUB_RHO_LIST, ax = ax1, dodge = True)

plt.savefig(outDir + outFileName, dpi = 400)

#Plot the true rho values on the y axis
for tick, text in zip(ax1.get_xticks(), ax1.get_xticklabels()):
    sample_name = text.get_text()  # "X" or "Y"

    #get the float value of rho corresponding with the tick
    rho_val = all_conf_ints_sub[all_conf_ints_sub['Sim_Rho'] == sample_name]
    rho_val = rho_val['Sim_float_rho'].unique()[0]

    # plot horizontal lines across the column, centered on the tick
    ax1.plot([tick-label_width/1.3, tick+label_width/1.3], [rho_val, rho_val],
            lw=1, color='r', linestyle = 'dotted')
    
# ax1.get_legend().remove()
# ax1.legend(bbox_to_anchor=(0.5, -0.205),
#             loc='upper center', ncol = 3, frameon = False,
#             columnspacing = 0.5, handletextpad = 0.1)
ax1.legend(bbox_to_anchor=(0.5, 0.99),
            loc='upper center', ncol = 3, frameon = False,
            columnspacing = 0, handletextpad = 0)
ax1.set_xticklabels(ax1.get_xticklabels(), rotation = 90)
ax1.set_xlabel(r'True Recombination Rate ($\rho$)')
ax1.set_title('Selection and Substructure')

###############################################################################
############################## Substructure Inset #############################
###############################################################################
# ax3 = fig.add_axes([0.55, 0.4, 0.08, 0.13])

# AF_THRESH_STR = ('20', '80')
# DIST_THRESH = 400
# fileTags = AF_THRESH_STR[0] + '_' + AF_THRESH_STR[1] + '_' +\
#           str(DIST_THRESH) + '.pkl'
# # all_D_vals_df = slim_util.make_dprime_df(dataDirSubInset, 
# #                     allele_freq_threshes = (0.2, 0.8), distance_thresh = 400,
# #                     four_hap_thresh = True, migration = True)

# all_D_vals_df = pd.read_pickle(dataDirSubInset + 'all_D_vals_df' + fileTags)
# all_D_vals_sel = pd.read_pickle(dataDirSelInset + 'all_D_vals_df' + fileTags)
# all_D_vals_sel['m_rate'] = 'Well Mixed'
# all_D_vals_df = pd.concat([all_D_vals_df, all_D_vals_sel], ignore_index = True)

# all_D_vals_df = slim_util.label_m_rate(all_D_vals_df)

# all_binned_rhos = []
# for name, group in all_D_vals_df.groupby(['m_rate', 'Rho']):
#     #Calculate the average D' value over distance
#     bin_ave_sim, binedges_sim, bin_nums_sim = binned_statistic(
#         group['Dist'].to_numpy(), 
#         group['d_prime'].to_numpy(), bins = 10)
#     rho_summaries = pd.DataFrame(zip(bin_ave_sim, binedges_sim[:-1]),
#                         columns = ['D_prime', 'Dist'])
#     rho_summaries['m_rate'] = name[0]
#     rho_summaries['Rho'] = name[1]
#     all_binned_rhos.append(rho_summaries)

# all_binned_rhos = pd.concat(all_binned_rhos, ignore_index = True)
# all_binned_rhos['m_rate'] = pd.Categorical(all_binned_rhos['m_rate'],
#                                 categories = ["Well Mixed", r"$m=10^{-2}$",
#                                 r"$m=10^{-3}$", r"$m=10^{-4}$", r"$m=10^{-5}$"])

# all_binned_rhos_04 = all_binned_rhos[all_binned_rhos['Rho'] == '1e-04']
# sns.lineplot(x = 'Dist', y = 'D_prime', data = all_binned_rhos_04, hue = 'm_rate',
#              style = 'Rho', palette = sns.color_palette("YlOrRd",
#                 n_colors= len(all_binned_rhos['m_rate'].unique())),
#                 linewidth = linewidth)

# ax3.yaxis.set_tick_params(pad = 0, length = 0.75)
# ax3.xaxis.set_tick_params(pad = 3, length = 0.75)
# ax3.set_ylabel('Linkage (D\')', labelpad=0.5)
# ax3.set_xlabel('Distance (bp)', labelpad=-1)
# ax3.set_ylim(0,1)
# ax3.get_legend().remove()

# all_binned_rhos_05 = all_binned_rhos[all_binned_rhos['Rho'] == '1e-05']
# ax4 = fig.add_axes([0.42, 0.71, 0.08, 0.13])
# sns.lineplot(x = 'Dist', y = 'D_prime', data = all_binned_rhos_05, hue = 'm_rate',
#              style = 'Rho', palette = sns.color_palette("YlOrRd",
#                 n_colors= len(all_binned_rhos['m_rate'].unique())),
#                 linewidth = linewidth, ax = ax4)

# ax4.yaxis.set_tick_params(pad = 0, length = 0.75)
# ax4.xaxis.set_tick_params(pad = 3, length = 0.75)
# ax4.set_ylabel('Linkage (D\')', labelpad=0.5)
# ax4.set_xlabel('Distance (bp)', labelpad=-1)
# ax4.set_ylim(0,1)
# ax4.get_legend().remove()


###############################################################################
########################## Varying NE Accuracy Panel ##########################
###############################################################################
vNE_RHO_LIST = [r"$2\times10^{-6}$",
                        r"$10^{-5}$",
                        r"$2\times10^{-5}$",
                        r"$10^{-4}$",
                        r"$2\times10^{-4}$",
                        r"$10^{-3}$"]
NUM_BOOTSTRAPS_vNE = 1000
NUM_GROUPS_vNE = 10

all_conf_ints_vNE = pd.read_pickle(dataDirVNE + "all_conf_ints_"+ \
                                        str(NUM_BOOTSTRAPS_vNE) + "_" + \
                                            str(NUM_GROUPS_vNE) + ".pkl")

#Relabel the Rho values and rename the accuracy column
all_conf_ints_vNE = slim_util.label_rho(all_conf_ints_vNE)
all_conf_ints_vNE['Estimated_Rho'] = all_conf_ints_vNE['est_rho']
all_conf_ints_vNE = all_conf_ints_vNE[
        all_conf_ints_vNE['Ne'].isin(['5000', '10000', '50000'])]
all_conf_ints_vNE = slim_util.label_Ne(all_conf_ints_vNE)
all_conf_ints_vNE['Ne'] = pd.Categorical(all_conf_ints_vNE['Ne'],
                                categories = [r"$N_e=5\times10^3$",
                                               r"$N_e=10^4$", 
                                               r"$N_e=5\times10^4$"])

print(all_conf_ints_vNE[all_conf_ints_vNE['Sim_Rho'] == r"$2\times10^{-6}$"])


#Now plot the accuracy
ax2 = axs[2]
sns.stripplot(x = 'Sim_Rho', y = 'Estimated_Rho',
                data = all_conf_ints_vNE, jitter = True, s = 2, alpha = 0.5,
                palette = sns.dark_palette("turquoise",
                n_colors= len(all_conf_ints_vNE['Ne'].unique())),
                hue = 'Ne', order = vNE_RHO_LIST, ax = ax2, dodge = True)

plt.savefig(outDir + outFileName, dpi = 400)

#Plot the true rho values on the y axis
for tick, text in zip(ax2.get_xticks(), ax2.get_xticklabels()):
    sample_name = text.get_text()  # "X" or "Y"

    #get the float value of rho corresponding with the tick
    rho_val = all_conf_ints_vNE[all_conf_ints_vNE['Sim_Rho'] == sample_name]
    rho_val = rho_val['Sim_float_rho'].unique()[0]

    # plot horizontal lines across the column, centered on the tick
    ax2.plot([tick-label_width/1.5, tick+label_width/1.5], [rho_val, rho_val],
            lw=linewidth, color='r', linestyle = 'dotted')

ax2.legend(frameon = False, columnspacing = 0, handletextpad = 0)
ax2.set_xticklabels(ax2.get_xticklabels(), rotation = 90)
ax2.set_xlabel(r'True Recombination Rate ($\rho$)')
ax2.set_title("Selection and Effective \n " + r'Population Size ($N_e$)')

plt.savefig(outDir + outFileName, dpi = 400)

