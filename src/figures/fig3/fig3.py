import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
#for running on desktop
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plot_neher as plne
import zaniniUtil as zu
from scipy import stats

NUM_BOOTSTRAPS = 1000

#For running on desktop
dataDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/"
vlDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/"
outDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/paper/fig3/"

#make the dataframe containg D' ratios
stat_df = zu.combine_drats(dataDir)
stat_df = zu.label_vl_drats(stat_df, vlDir)

#filter the d' ratio dataframe
stat_df = stat_df[stat_df['Fragment'] != 'F5']
stat_df = stat_df[stat_df['Time_Diff'].gt(50)]
stat_df = stat_df[stat_df['Participant'] != 'p10']


all_par_ests = []
all_par_fits = []
x_vals = stat_df['Dist_X_Time'].unique()

#Estimate Rates Specifically for each individual
for name, group in stat_df.groupby('Participant'):
    print(group)
    print(name)
    #get estimates for the individual
    lower_fit, upper_fit, estimate_df = plne.bootstrap_rho(group, NUM_BOOTSTRAPS)
    estimate_df['Participant'] = name
    binned_rat, binedges, bin_nums = stats.binned_statistic(
        group['Dist_X_Time'].to_numpy(), 
        group['d_ratio'].to_numpy(), bins = 100)


    mid_fit = np.quantile(estimate_df['Estimated_Rho'], 0.5)
    mid_fit = estimate_df.iloc[(estimate_df['Estimated_Rho']-mid_fit).abs().argsort()[:2]]

    fit_vals_low = [plne.neher_leitner(x, lower_fit['c0'].to_numpy()[0], lower_fit['c1'].to_numpy()[0], lower_fit['c2'].to_numpy()[0]) for x in binedges]
    fit_vals_high = [plne.neher_leitner(x, upper_fit['c0'].to_numpy()[0], upper_fit['c1'].to_numpy()[0], upper_fit['c2'].to_numpy()[0]) for x in binedges]
    fit_vals_mid = [plne.neher_leitner(x, mid_fit['c0'].to_numpy()[0], mid_fit['c1'].to_numpy()[0], mid_fit['c2'].to_numpy()[0]) for x in binedges]

    par_fits = pd.DataFrame(zip(binned_rat, binedges, fit_vals_low, fit_vals_mid, fit_vals_high),
                    columns=['ratio_bins', 'bin_edges', 'lower_conf', 'mid_conf', 'high_conf'])
    par_fits['Participant'] = name

    all_par_fits.append(par_fits)
    all_par_ests.append(estimate_df)

    
all_par_ests = pd.concat(all_par_ests, ignore_index=True)
all_par_fits = pd.concat(all_par_fits, ignore_index=True)
print(all_par_fits)


################################ Estimate Plot ################################
plt.rcParams.update({
  "text.usetex": True
})
sns.set(rc={'figure.figsize':(20,7.5)}, font_scale = 2)
sns.set_palette("tab10")
ax = sns.boxplot(x = 'Participant', y = 'Estimated_Rho', data = all_par_ests,
    order = ['p1', 'p2', 'p3', 'p4', 'p5', 'p7', 'p8', 'p9', 'p11'])
ax.axhline(0.000008, linestyle = 'dashed', color = 'tab:green')
ax.axhline(0.000014, color = 'tab:green')
ax.axhline(0.00002, linestyle = 'dashed', color = 'tab:green')
ax.set_ylim(0, 0.0002)


plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
plt.xlabel(r'Participant')
plt.ylabel(r'Estimated Value of $\rho$')

plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
plt.tight_layout()
plt.savefig(outDir + "fig3.jpg")
plt.close()



