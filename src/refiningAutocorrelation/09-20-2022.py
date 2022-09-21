from statistics import median
import sys
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import plot_neher as plne
from scipy import optimize

#Today I am trying to show why we are underestimating low recombination rates
outdir = outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/methodValidation/09-20-2022/'
NOISE_LEVEL = 0.01

# ############################## Scenario 1 #####################################
#high rho
all_estimate_df = []
for i in range(100):
    print(i)
    X_MAX = 100000
    STEP = 10000
    SUBSET_THRESHOLDS = np.flip(np.arange(10000, X_MAX + STEP, STEP))
    C0 = 0.716570
    C1 = 0.518331
    C2 = 0.000255
    NUM_POINTS = 200


    #The functional form we are fitting to the data
    def known_equation(x, c0, c1, c2):
        return c0 + c1 * (1-np.exp(-c2 * x))

    #Get x values for samples
    mydata = np.random.uniform(0, X_MAX, NUM_POINTS)


    #Make data that follows the equation but add random noise
    mapped_data = np.apply_along_axis(known_equation, 0, mydata, c0=C0 , c1 = C1, c2 =C2)
    noise = np.random.normal(0, NOISE_LEVEL, NUM_POINTS)
    mapped_data = mapped_data + noise

    #Put the results into a neat dataframe and plot them
    result_df = pd.DataFrame(list(zip(mydata, mapped_data)), columns = ['x', 'y'])
    sns.scatterplot(x = 'x', y = 'y', data = result_df)


    all_fits_df = []
    estimate_df = []

    #Now fit to varying subsets of the data
    for curr_thresh in SUBSET_THRESHOLDS:
        #Get the data below the threshold
        curr_data = result_df[result_df['x'] < curr_thresh]
        if len(curr_data) < 3:
            continue

        #Plot so we can see what data ends up in each group
        sns.scatterplot(x = 'x', y = 'y', data = curr_data, label = curr_thresh)

        #Fit the curve
        coeffs, fit_dat = optimize.curve_fit(plne.neher_leitner, 
        curr_data['x'], curr_data['y'],p0 = [0, 0.26, .0000439], 
        maxfev = 10000)

        fit_vals = [plne.neher_leitner(x, coeffs[0], coeffs[1], coeffs[2])
                    for x in np.arange(0, X_MAX, STEP)]
        fit_vals = pd.DataFrame(list(zip(fit_vals, np.arange(0, X_MAX, STEP))), columns = ['y', 'x'])
        fit_vals['subset'] = curr_thresh
        
        estimate_df.append([curr_thresh, coeffs[0], coeffs[1], coeffs[2], coeffs[1] * coeffs[2]])
        all_fits_df.append(fit_vals)

    all_fits_df = pd.concat(all_fits_df, ignore_index=True)
    estimate_df = pd.DataFrame(estimate_df, columns = ['subset', 'c0', 'c1', 'c2', 'd\'_0'])
    

    plt.savefig(outdir + 'plotting_subsets_1_reduced_noise.png')
    plt.close()

    sns.lineplot(x = 'x', y = 'y', data = all_fits_df, hue = 'subset')
    sns.scatterplot(x = 'x', y = 'y', data = result_df)
    plt.savefig(outdir + 'subset_fits_1_reduced_noise.png')
    plt.close()

    estimate_df['rep'] = i
    all_estimate_df.append(estimate_df)


all_estimate_df = pd.concat(all_estimate_df, ignore_index=True)
print(all_estimate_df)
ax = sns.stripplot(all_estimate_df['subset'], all_estimate_df['d\'_0'])

# distance across the "X" or "Y" stipplot column to span, in this case 40%
median_width = 0.4

for tick, text in zip(ax.get_xticks(), ax.get_xticklabels()):
    sample_name = text.get_text()  # "X" or "Y"
    # calculate the median value for all replicates of either X or Y
    subset_data = all_estimate_df[all_estimate_df['subset']==float(sample_name)]
    median_val = np.median(subset_data['d\'_0'])


    # plot horizontal lines across the column, centered on the tick
    ax.plot([tick-median_width/2, tick+median_width/2], [median_val, median_val],
            lw=4, color='k')
plt.axhline(y=C1 * C2, color='r')
plt.ylim(0, 0.001)
plt.savefig(outdir + 'stripplot_1_reduced_noise.png')
plt.close()

############################## Scenario 2 #####################################
#high rho
all_estimate_df = []
for i in range(100):
    print(i)
    X_MAX = 100000
    STEP = 10000
    SUBSET_THRESHOLDS = np.flip(np.arange(10000, X_MAX + STEP, STEP))
    C0 = 0.119403
    C1 = 1.608203
    C2 = 1.288099e-06
    NUM_POINTS = 200


    #The functional form we are fitting to the data
    def known_equation(x, c0, c1, c2):
        return c0 + c1 * (1-np.exp(-c2 * x))

    #Get x values for samples
    mydata = np.random.uniform(0, X_MAX, NUM_POINTS)


    #Make data that follows the equation but add random noise
    mapped_data = np.apply_along_axis(known_equation, 0, mydata, c0=C0 , c1 = C1, c2 =C2)
    noise = np.random.normal(0, NOISE_LEVEL, NUM_POINTS)
    mapped_data = mapped_data + noise

    #Put the results into a neat dataframe and plot them
    result_df = pd.DataFrame(list(zip(mydata, mapped_data)), columns = ['x', 'y'])
    sns.scatterplot(x = 'x', y = 'y', data = result_df)


    all_fits_df = []
    estimate_df = []

    #Now fit to varying subsets of the data
    for curr_thresh in SUBSET_THRESHOLDS:
        #Get the data below the threshold
        curr_data = result_df[result_df['x'] < curr_thresh]
        if len(curr_data) < 3:
            continue

        #Plot so we can see what data ends up in each group
        sns.scatterplot(x = 'x', y = 'y', data = curr_data, label = curr_thresh)

        #Fit the curve
        coeffs, fit_dat = optimize.curve_fit(plne.neher_leitner, 
        curr_data['x'], curr_data['y'],p0 = [0, 0.26, .0000439], 
        maxfev = 10000)

        fit_vals = [plne.neher_leitner(x, coeffs[0], coeffs[1], coeffs[2])
                    for x in np.arange(0, X_MAX + STEP, STEP)]
        fit_vals = pd.DataFrame(list(zip(fit_vals, np.arange(0, X_MAX + STEP, STEP))), columns = ['y', 'x'])
        fit_vals['subset'] = curr_thresh
        
        estimate_df.append([curr_thresh, coeffs[0], coeffs[1], coeffs[2], coeffs[1] * coeffs[2]])
        all_fits_df.append(fit_vals)

    all_fits_df = pd.concat(all_fits_df, ignore_index=True)
    estimate_df = pd.DataFrame(estimate_df, columns = ['subset', 'c0', 'c1', 'c2', 'd\'_0'])
    

    plt.savefig(outdir + 'plotting_subsets_2_reduced_noise.png')
    plt.close()

    sns.lineplot(x = 'x', y = 'y', data = all_fits_df, hue = 'subset')
    sns.scatterplot(x = 'x', y = 'y', data = result_df)
    plt.savefig(outdir + 'subset_fits_2_reduced_noise.png')
    plt.close()

    estimate_df['rep'] = i
    all_estimate_df.append(estimate_df)


all_estimate_df = pd.concat(all_estimate_df, ignore_index=True)
ax = sns.stripplot(all_estimate_df['subset'], all_estimate_df['d\'_0'])

# distance across the "X" or "Y" stipplot column to span, in this case 40%
median_width = 0.4

for tick, text in zip(ax.get_xticks(), ax.get_xticklabels()):
    sample_name = text.get_text()  # "X" or "Y"
    # calculate the median value for all replicates of either X or Y
    subset_data = all_estimate_df[all_estimate_df['subset']==float(sample_name)]
    median_val = np.median(subset_data['d\'_0'])

    # plot horizontal lines across the column, centered on the tick
    ax.plot([tick-median_width/2, tick+median_width/2], [median_val, median_val],
            lw=4, color='k')
plt.axhline(y=C1 * C2, color='r')
plt.ylim(0, 0.00001)
plt.savefig(outdir + 'stripplot_2_reduced_noise.png')
plt.close()

############################## Scenario 2 #####################################
#high rho with half the data 
all_estimate_df = []
for i in range(100):
    print(i)
    X_MAX = 50000
    STEP = 10000
    SUBSET_THRESHOLDS = np.flip(np.arange(10000, X_MAX + STEP, STEP))
    C0 = 0.119403
    C1 = 1.608203
    C2 = 1.288099e-06
    NUM_POINTS = 200


    #The functional form we are fitting to the data
    def known_equation(x, c0, c1, c2):
        return c0 + c1 * (1-np.exp(-c2 * x))

    #Get x values for samples
    mydata = np.random.uniform(0, X_MAX, NUM_POINTS)


    #Make data that follows the equation but add random noise
    mapped_data = np.apply_along_axis(known_equation, 0, mydata, c0=C0 , c1 = C1, c2 =C2)
    noise = np.random.normal(0, NOISE_LEVEL, NUM_POINTS)
    mapped_data = mapped_data + noise

    #Put the results into a neat dataframe and plot them
    result_df = pd.DataFrame(list(zip(mydata, mapped_data)), columns = ['x', 'y'])
    sns.scatterplot(x = 'x', y = 'y', data = result_df)


    all_fits_df = []
    estimate_df = []

    #Now fit to varying subsets of the data
    for curr_thresh in SUBSET_THRESHOLDS:
        #Get the data below the threshold
        curr_data = result_df[result_df['x'] < curr_thresh]
        if len(curr_data) < 3:
            continue

        #Plot so we can see what data ends up in each group
        sns.scatterplot(x = 'x', y = 'y', data = curr_data, label = curr_thresh)

        #Fit the curve
        coeffs, fit_dat = optimize.curve_fit(plne.neher_leitner, 
        curr_data['x'], curr_data['y'],p0 = [0, 0.26, .0000439], 
        maxfev = 10000)

        fit_vals = [plne.neher_leitner(x, coeffs[0], coeffs[1], coeffs[2])
                    for x in np.arange(0, X_MAX + STEP, STEP)]
        fit_vals = pd.DataFrame(list(zip(fit_vals, np.arange(0, X_MAX + STEP, STEP))), columns = ['y', 'x'])
        fit_vals['subset'] = curr_thresh
        
        estimate_df.append([curr_thresh, coeffs[0], coeffs[1], coeffs[2], coeffs[1] * coeffs[2]])
        all_fits_df.append(fit_vals)

    all_fits_df = pd.concat(all_fits_df, ignore_index=True)
    estimate_df = pd.DataFrame(estimate_df, columns = ['subset', 'c0', 'c1', 'c2', 'd\'_0'])
    

    plt.savefig(outdir + 'plotting_subsets_2_50Samples_reduced_noise.png')
    plt.close()

    sns.lineplot(x = 'x', y = 'y', data = all_fits_df, hue = 'subset')
    sns.scatterplot(x = 'x', y = 'y', data = result_df)
    plt.savefig(outdir + 'subset_fits_2_50Samples_reduced_noise.png')
    plt.close()

    estimate_df['rep'] = i
    all_estimate_df.append(estimate_df)


all_estimate_df = pd.concat(all_estimate_df, ignore_index=True)
ax = sns.stripplot(all_estimate_df['subset'], all_estimate_df['d\'_0'])

# distance across the "X" or "Y" stipplot column to span, in this case 40%
median_width = 0.4

for tick, text in zip(ax.get_xticks(), ax.get_xticklabels()):
    sample_name = text.get_text()  # "X" or "Y"
    # calculate the median value for all replicates of either X or Y
    subset_data = all_estimate_df[all_estimate_df['subset']==float(sample_name)]
    median_val = np.median(subset_data['d\'_0'])

    # plot horizontal lines across the column, centered on the tick
    ax.plot([tick-median_width/2, tick+median_width/2], [median_val, median_val],
            lw=4, color='k')
plt.axhline(y=C1 * C2, color='r')
plt.ylim(0, 0.00001)
plt.savefig(outdir + 'stripplot_2_50Samples_reduced_noise.png')
plt.close()