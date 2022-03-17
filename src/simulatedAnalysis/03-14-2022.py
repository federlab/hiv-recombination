import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
# #For running on Desktop
# sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import os
import numpy as np
import pandas as pd
import seaborn as sns
import plot_neher as plne
import matplotlib.pyplot as plt
import zaniniUtil as zu

#I am going to use this script to make figures for the Harris Lab Meeting

#For running on cluster
dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_02_24/'
outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_02_24/'

# #For running on desktop
# dataDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/slimDatasets/2022_02_24/'
# outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/slimDatasets/2022_02_24/'

if not os.path.isdir(outDir):
    os.mkdir(outDir)

###############################################################################
#Make dataframes to store the resulting frequencies in
all_dt_freqs = []
all_d_freqs = []
estimate_df = []

for curr_dataset in os.listdir(dataDir):
    if curr_dataset.split('_')[0] != 'mu1e-05':
        continue
    curr_rho = curr_dataset.split('_')[1]
    curr_rho = curr_rho.split('o')[1]

    #Get the current data
    recombination_df = pd.read_pickle(
                        dataDir + curr_dataset + "/neher_res/recombination")
    recombination_df['dist'] = recombination_df['Locus_2'] - \
                                recombination_df['Locus_1']
    recombination_df['Dist_x_Time'] = (recombination_df['Curr_Timepoint'] - \
                recombination_df['Last_Timepoint']) * recombination_df['dist']

    mutation_df = pd.read_pickle(
                        dataDir + curr_dataset + "/neher_res/mutation")
    mutation_df['Dist_x_Time'] = mutation_df['Dist_x_Time'].abs()

    ###### Plot by Dist_x_Time ######
    #make the set of bins
    dtbinset = [(x, x + 2500) for x in range(0, 55000, 2500)]
    dt_freqs = plne.bin_curve(
        recombination_df, mutation_df, dtbinset, 'Dist_x_Time')
    dt_freqs['Dataset'] = curr_dataset
    dt_freqs['Rho'] = curr_rho
    all_dt_freqs.append(dt_freqs)

    ###### Plot by dist ######    
    #make the set of bins
    dbinset = [(x, x + 50) for x in range(0, 500, 50)]
    d_freqs = plne.bin_curve(
        recombination_df, mutation_df, dbinset, 'dist')
    d_freqs['Dataset'] = curr_dataset
    d_freqs['Rho'] = curr_rho
    all_d_freqs.append(d_freqs)

    # ########################### Fitting #######################################
    #perform our fitting
    try:
        coeffs, fit_df = plne.run_neher_fit(c0_fixed = False, lower_bounds = [0,0,0],
                                    upper_bounds = [1,1,1], 
                                    initial = [0.13, 0.1804, 0.000114],
                                    test_results = dt_freqs)
    except ValueError as my_err:
        print(my_err, file = sys.stderr)
        coeffs = [0,0,0]
        fit_df = pd.DataFrame(columns= ['x_vals', 'fitted_vals'], index = range(1))

    curr_estimate = plne.estimate_recombination_rate(c0 = coeffs[0], c1 = coeffs[1], c2 = coeffs[2])
    estimate_df.append([curr_estimate, curr_dataset, curr_rho])

#Concatenate all of the frequencies together
all_d_freqs = pd.concat(all_d_freqs, ignore_index = True )
all_dt_freqs = pd.concat(all_dt_freqs, ignore_index= True)

estimate_df = pd.DataFrame(estimate_df, columns=["Est_Rho", 'Dataset', 'Sim_Rho'])

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

#add jitter to the rho values
x_vals = np.linspace(0.00001, 0.000105, 10)
def jitter(values,j):
    return values + np.random.normal(j,0.1,values.shape)
estimate_df['Sim_int_rho']

print(estimate_df)
sns.set(rc={'figure.figsize':(10,5)}, font_scale = 1)
print(estimate_df['Sim_int_rho'])


ax = sns.stripplot(x = 'Sim_Rho', y = 'Est_Rho', data = estimate_df, jitter = True,
    order = ["0.001", "2e-04", "1e-04", "2e-05", "1e-05", "2e-06"])

# distance across the "X" or "Y" stipplot column to span, in this case 40%
label_width = 0.4

for tick, text in zip(ax.get_xticks(), ax.get_xticklabels()):
    sample_name = text.get_text()  # "X" or "Y"

    # calculate the median value for all replicates of either X or Y
    rho_val = rhoDict[sample_name]

    # plot horizontal lines across the column, centered on the tick
    ax.plot([tick-label_width/2, tick+label_width/2], [rho_val, rho_val],
            lw=2, color='k')

plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
plt.xlabel("Simulation Value of Rho")
plt.ylabel("Estimated Value of Rho")
plt.ylim(0.0000001, 0.1)
plt.yscale('log')
plt.tight_layout()
plt.savefig(outDir + "comparedEstimates_stripplot_all_rhos.jpg")
plt.close()


