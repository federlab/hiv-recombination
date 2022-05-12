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

#This script just makes a basic plot of autoccorelation estimates for simulated data against
#the actual values in the simulation
#from this analysis, it looks like fitting on the sampled data doesn't quite work
THRESHOLD = 0.2
DIST_TIME_MAX = 50000
NUM_BINS = 25
NUM_PER_BIN = 500 #number of observations per bin (sampled with replacement)

#make a set of bins and sample similar amounts in every bin
bin_starts = range(0, DIST_TIME_MAX, int(DIST_TIME_MAX/NUM_BINS))
bin_edges = [(x, x + int(DIST_TIME_MAX/NUM_BINS)) for x in bin_starts]

#For running on cluster
dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_04_20/'
outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_04_20/'

dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_04_20/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_04_20/05-04-2022/'


# short_times_all = []
# long_times_all = []
# all_first_bin = []
# all_middle_bin = []

# for curr_data in os.listdir(dataDir):
#     print(curr_data)

#     #only get the data directories, not hidden files
#     if curr_data[0] == '.':
#         continue
#     run_info = curr_data.split('_')
#     sim_rho = run_info[1]
#     sim_rho = sim_rho[3:]
#     rep = run_info[-1]

#     linkage_file = dataDir + curr_data + "/linkage/r2_and_D"
#     d_ratio_file = dataDir + curr_data + "/linkage/d_ratio"

#     if os.path.exists(d_ratio_file):
#         stat_df = pd.read_pickle(d_ratio_file)
#     else:
#         stat_df = autocorr.calculate_d_ratios(linkage_file, THRESHOLD)
#         stat_df.to_pickle(d_ratio_file)

#     stat_df['dist'] = stat_df['Dist_X_Time']/stat_df['Time_Diff']

#     short_times  = stat_df[stat_df['Time_Diff'] < 10]
#     long_times = stat_df[stat_df['Time_Diff'] > 10]
#     short_times_all.append(short_times)
#     long_times_all.append(long_times)

#     first_bin = stat_df[stat_df['Dist_X_Time'].between(0, 2500)]
#     all_first_bin.append(first_bin)
#     middle_bin = stat_df[stat_df['Dist_X_Time'].between(7500, 10000)]
#     all_middle_bin.append(middle_bin)

# short_times_all = pd.concat(short_times_all, ignore_index= True)
# long_times_all = pd.concat(long_times_all, ignore_index= True)
# all_first_bin = pd.concat(all_first_bin, ignore_index= True)
# all_middle_bin = pd.concat(all_middle_bin, ignore_index= True)


# #plot all of the D values
# sns.scatterplot(x = 'Dist_X_Time', y = 'd_i', data = short_times_all, alpha = 0.005)
# plt.savefig(outDir + "d_i_short_times.jpg")
# plt.close()

# sns.scatterplot(x = 'Dist_X_Time', y = 'd_i', data = long_times_all, alpha = 0.005)
# plt.xlim(0, 2500)
# plt.savefig(outDir + "d_i_long_times.jpg")
# plt.close()

# sns.scatterplot(x = 'Dist_X_Time', y = 'd_i', data = all_first_bin, alpha = 0.005)
# plt.savefig(outDir + "d_i_first_bin.jpg")
# plt.close()

# sns.scatterplot(x = 'Dist_X_Time', y = 'd_i', data = all_middle_bin, alpha = 0.005)
# plt.savefig(outDir + "d_i_middle_bin.jpg")
# plt.close()

possible_D_i = [x * 0.1 for x in range(2,11)]
possible_D_i_1 = [x * 0.1 for x in range(2,11)]

my_results = []

for j in range(len(possible_D_i)):
    for k in range(len(possible_D_i_1)):
        d_i = possible_D_i[j]
        d_i_1 = possible_D_i_1[k]

        result = -np.log(d_i_1/d_i)
        if d_i > d_i_1:
            my_results.append([result, d_i, d_i_1, 'decrease'])
        else: my_results.append([result, d_i, d_i_1, 'increase'])

my_results = pd.DataFrame(my_results, columns= ['d_rat', 'd_i', 'd_i_1', 'change'])


sns.scatterplot(x = 'd_i', y = 'd_rat', data = my_results, hue = 'change')
plt.savefig(outDir + "possible_d_rat.jpg")
plt.close()