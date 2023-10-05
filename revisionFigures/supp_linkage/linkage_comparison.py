import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import pandas as pd
import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import slimUtil as slim_util
import zaniniUtil as zanini_util
from scipy.stats import binned_statistic
from matplotlib import rcParams

#This file will make a supplementary figure that shows the linkage (D' vs distance) for the in vivo data
# and the simulation data. 

#In this file I will be checking the linkage at long range and comparing it to p4 and p7
dataDirSub = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/revision/simulated_linkage_sub/'
dataDirSel = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/revision/simulated_linkage_sel/'
inVivoDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/revision/supp_figs/'

###############################################################################
##################### Setting the figure parameters ###########################
###############################################################################
#plot the estimates to show how accurate they are
params = {'figure.figsize':(7, 3), 'axes.labelsize': 6,'axes.titlesize':6,  
          'legend.fontsize': 6, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
          'legend.title_fontsize': 6}
rcParams.update(params)

fig, axs = plt.subplots(1, 2, sharey= True)
plt.subplots_adjust(wspace = 0.05, hspace=0)

linewidth = 1

###############################################################################
########################## Simulated Processing ###############################
###############################################################################
LOW_FREQ = '20'
HIGH_FREQ = '80'
DIST_THRESH = '400'

#Read in the substructured data
fileTags = LOW_FREQ + '_' + HIGH_FREQ + '_' + DIST_THRESH + '.pkl'
all_D_vals_df = pd.read_pickle(dataDirSub + 'all_D_vals_df' + fileTags)

#Read in the well mixed data
all_D_vals_df_sel = pd.read_pickle(dataDirSel + 'all_D_vals_df' + fileTags)
all_D_vals_df_sel['m_rate'] = 'Well Mixed'

all_D_vals_df  = pd.concat(
                        [all_D_vals_df , all_D_vals_df_sel], ignore_index=True)
all_D_vals_df  = slim_util.label_m_rate(all_D_vals_df )

all_binned_rhos_sim = []
for name, group in all_D_vals_df.groupby(['m_rate', 'Rho']):
    #Calculate the average D' value over distance
    bin_ave_sim, binedges_sim, bin_nums_sim = binned_statistic(
        group['Dist'].to_numpy(), 
        group['d_prime'].to_numpy(), bins = 10)
    rho_summaries = pd.DataFrame(zip(bin_ave_sim, binedges_sim[:-1]),
                        columns = ['D_prime', 'Dist'])
    rho_summaries['m_rate'] = name[0]
    rho_summaries['Sim_Rho'] = name[1]
    all_binned_rhos_sim.append(rho_summaries)


all_binned_rhos_sim= pd.concat(all_binned_rhos_sim, ignore_index = True)
all_binned_rhos_sim['m_rate'] = pd.Categorical(all_binned_rhos_sim['m_rate'],
                                categories = ["Well Mixed", r"$m=10^{-2}$",
                                r"$m=10^{-3}$", r"$m=10^{-4}$", r"$m=10^{-5}$"])
all_binned_rhos_sim = slim_util.label_rho(all_binned_rhos_sim)

###############################################################################
############################ In Vivo Processing ###############################
###############################################################################

#We need to call the in vivo version of the function
in_vivo_D = zanini_util.make_dprime_df(inVivoDir, 
                        allele_freq_threshes = (0.2, 0.8), 
                        distance_thresh = 400, four_hap_thresh = True,
                        frags_to_exclude = ['F5'])

label_Dict = {
        'p4': 'Superinfection',
        'p7': 'Superinfection',
        'p3': 'Multi Founder',
        'p10': 'Multi Founder',
        'p11': 'Single Founder',
        'p9': 'Single Founder',
        'p8': 'Single Founder',
        'p6': 'Single Founder',
        'p5': 'Single Founder',
        'p2': 'Single Founder',
        'p1': 'Single Founder'}

in_vivo_D['inf_type'] = in_vivo_D['Participant'].map(label_Dict)

all_binned_rhos_inv = []
#Calculate the average D' value over distance
for name, group in in_vivo_D.groupby(['inf_type']):
    #Calculate the average D' value over distance
    bin_ave_sim, binedges_sim, bin_nums_sim = binned_statistic(
        np.array(group['Dist'].tolist()), 
        np.array(group['d_prime'].tolist()), bins = 10)
    rho_summaries = pd.DataFrame(zip(bin_ave_sim, binedges_sim[:-1]),
                        columns = ['D_prime', 'Dist'])
    rho_summaries['inf_type'] = name
    all_binned_rhos_inv.append(rho_summaries)

all_binned_rhos_inv= pd.concat(all_binned_rhos_inv, ignore_index = True)

###############################################################################
############################# Plotting the Data ###############################
###############################################################################
iter_ax = axs.flatten()

for i in range(len(all_binned_rhos_sim['Sim_Rho'].unique())):
    curr_ax = iter_ax[i]
    curr_rho = all_binned_rhos_sim['Sim_Rho'].unique()[i]
    curr_binned_rhos_sim = all_binned_rhos_sim[
                            all_binned_rhos_sim['Sim_Rho'] == curr_rho]

    sns.lineplot(x='Dist', y='D_prime', hue='m_rate', 
                data=curr_binned_rhos_sim , 
                palette = sns.color_palette("YlOrRd", 
                n_colors= len(curr_binned_rhos_sim ['m_rate'].unique())),
                linewidth = linewidth, ax = curr_ax)

    sns.lineplot(x='Dist', y='D_prime', hue='inf_type', 
                data=all_binned_rhos_inv, palette = sns.color_palette("BuPu",
                n_colors= len(all_binned_rhos_inv['inf_type'].unique())),
                linewidth = linewidth, ax = curr_ax, linestyle = "--")
    
    curr_ax.set_ylim(0, 1)
    curr_ax.set_xlabel('Distance (bp)')
    curr_ax.set_title(r'Simulated $\rho$ = ' + str(curr_rho))

handles, labels = axs[1].get_legend_handles_labels()


#Separate out and make the legend for the simulated data
sim_handles = handles[:5]
sim_labels = labels[:5]
axs[1].legend(handles = sim_handles, labels = sim_labels, title = r'Simulated Migration Rate', ncol = 2, columnspacing = 1, handlelength = 1.5)

#Separate out and make the legend for the in vivo data
inv_handles = handles[5:]
inv_labels = labels[5:]
axs[0].legend(handles = inv_handles, labels = inv_labels, title = r'$\it{In\ Vivo}$ Characteristics', ncol = 3, columnspacing = 1, handlelength = 1.5)


axs[0].set_ylabel("Linkage (Average D\')")
plt.savefig(outDir + "linkage.jpg", dpi = 400)