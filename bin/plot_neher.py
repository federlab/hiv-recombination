import sys
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import optimize

#This script is specifically for plotting the results of the neher and
#leitner analysis

############ Helper Functions for plotting and estimating the rate ############
def bin_curve(recombination_df, mutation_df, bins, bin_type):
    """ This script calculates the frequencies of mutation and recombination
    test successes in each given bin.
    ---------------------------------------------------------------------------
    Params
    ------------
    recombination_df : pd.dataframe, one column indicates whether a pair of 
                loci triggered a three haplotype test. Two columns give the
                loci pair and a third column indicates whether they passed the  
                test. Additionally it has two columns indicating what 
                timepoints these were observed between. The final column 
                indicates the support of the test. 
    mutation_df : pd.dataframe, dataframe similar to recombination df. 
                Indicates when a test has been conducted to see if a new 
                mutation has been observed. Also indicates the results and 
                support of the test.
    bins : list of tuples, a list of tuples where the first element in the 
                tuple is the bin start and the second element in the tuple is
                the bin end. Tests are grouped in these bins to calculate the
                frequency of success to be used for plotting and curve fitting 
    bin_type : string, the column name of the column used for the binning. This
                is generally 'dist' or 'Dist_x_Time'

    Returns
    -------------
    dataset_freqs : pd.dataframe, a dataframe with the 'window' column which
                indicates what bin the tests are in and then additionally 
                columns listing the number of tests and the number of successes
                for both the recombination and mutation tests. Additionally,
                there are two error columns which are just the square root of
                the number of tests.
    """
    dataset_freqs = []

    for currBin in bins:
        bin_start = currBin[0]
        bin_end = currBin[1]
        #added 1/14/22 plot the bins at the center instead of the start
        bin_center = int((currBin[1] - currBin[0])/2) + bin_start

        #make a place to store the frequencies of observing events
        mut_frequencies = []
        mut_successes = []
        recomb_frequencies = []
        recomb_successes = []
        mut_tests = []
        recomb_tests = []

        #get all of the datapoints in our bin
        curr_recomb = recombination_df[recombination_df[bin_type].between(
                bin_start, bin_end)]
        curr_mut = mutation_df[mutation_df[bin_type].between(
                bin_start, bin_end)]


        #Calculate the frequencies in each bin
        if curr_recomb.shape[0] > 0:
            recomb_true = curr_recomb[curr_recomb['Test_Passed'] == True]
            recomb_frequencies.append(recomb_true.shape[0]/curr_recomb.shape[0])
            recomb_successes.append(recomb_true.shape[0])
            recomb_tests.append(curr_recomb.shape[0])

        else: 
            recomb_frequencies.append(0)
            recomb_successes.append(0)
            recomb_tests.append(0)

        if curr_mut.shape[0] > 0:
            mut_true = curr_mut[curr_mut['Test_Passed'] == True]
            mut_frequencies.append(mut_true.shape[0]/curr_mut.shape[0])
            mut_successes.append(mut_true.shape[0])
            mut_tests.append(curr_mut.shape[0])
        else: 
            mut_frequencies.append(0)
            mut_successes.append(0)
            mut_tests.append(0)

        curr_frequencies = pd.DataFrame(list(zip(mut_frequencies,recomb_frequencies)),
                            columns = ['mut_frequencies', 'recomb_frequencies'])
        curr_frequencies['window'] = bin_end
        curr_frequencies['Mutation Tests'] = mut_tests
        curr_frequencies['Recombination Tests'] = recomb_tests
        curr_frequencies['Mutation Successes'] = mut_successes
        curr_frequencies['Recombination Successes'] = recomb_successes
        dataset_freqs.append(curr_frequencies)

    dataset_freqs = pd.concat(dataset_freqs, ignore_index= True)

    #get the error bars for each set of tests
    dataset_freqs['Mut Error'] = 1 / np.sqrt(dataset_freqs['Mutation Tests'])
    dataset_freqs['Recomb Error'] = 1 / np.sqrt(dataset_freqs['Recombination Tests'])

    return dataset_freqs

def estimate_recombination_rate(c0, c1, c2, empirical = False, MpApB = 0):
    """ Takes in the coefficients from fitting the curve to our data. Then
    returns our estimate of the recombination rate
    ---------------------------------------------------------------------------
    Params
    ------------
    c0 : float, the intercept of our curve fit
    c1 : float, coefficient for second term
    c2 : float, coefficient in the exp
    empirical : bool, whether to use the empirical MpApB
    MpApB : float, value of empirical MpApB to use

    Returns
    -------------
    rec_rate : the per virus recombination rate
    """
    numerator = c1 * c2 
    if empirical:
        if MpApB == 0:
            raise ValueError("Must specify MpApB")
        denominator = MpApB - c0  
    else:
        denominator = np.log(1/(1-c0))
        denominator = np.log(1/(1 - c0 - c1)) - denominator
        denominator = (1-c0) * denominator
    return numerator/ denominator

######################### Functionality for Fitting ###########################
def neher_leitner(dDeltaT, c0, c1, c2):
    """The function Neher and Leitner fit to determine recombination rate"""
    return (c0 + (c1 * (1 - np.exp(-c2 * dDeltaT))))

def neher_leitner_c0_fixed(dDeltaT, c1, c2):
    """The function Neher and Leitner fit to determine recombination rate"""
    return ((c1 * (1 - np.exp(-c2 * dDeltaT))))

def residual(x0, dDeltaT, data, rec_error, fixed_c0 = None):
    """ Calculate the residuals for our fit
    ---------------------------------------------------------------------------
    Params
    ------------
    x0 :      list, initial guesses at parameters c1 and c2. Also contains c0
              if it is not fixed
    dDeltaT : pd.df column, the 'window' column of our dataframe
    data :    pd.df column, the column our our data with the frequency of test
              success for our recombination tests
    rec_error:pd.df column, the column with the weights for the fit (based
              on the number of tests run)
    fixed_c0: float, the value of c0 if c0 is fixed rather than being fitted.
              If c0 is fitted, then this should be None.

    Returns
    -------------
    weighted_resids : np.array the residuals taking into account the model
    """
    #set the values depending on whether c0 is fixed or not
    if fixed_c0 is not None and len(x0) == 2:
        c0 = fixed_c0
        c1 = x0[0]
        c2 = x0[1]
    elif len(x0) == 3:
        c0 = x0[0]
        c1 = x0[1]
        c2 = x0[2]
    else:
        raise ValueError('Specification of c0, c1, c2 was not correct.')

    #this is the equation we are using for our fit
    model = neher_leitner(dDeltaT, c0, c1, c2)
    resids = data - model
    weighted_resids = resids * (1 + rec_error)
    return weighted_resids

def run_neher_fit(c0_fixed, lower_bounds, upper_bounds, initial, test_results):
    """ Performs fitting of data to a curve using the functional form described
    by Neher and Leitner. Returns one list containing the estimates for 
    the coefficients and a dataframe with fitted values for plotting.
    ---------------------------------------------------------------------------
    Params
    ------------
    c0_fixed :      bool, whether the constant c0 will be fixed or fitted
    lower_bounds :  list, A list of floats containing the minimum values to try
                    for c0, c1 and c2. (c0 not included if it is fixed)
    upper_bounds :  list, A list of floats containing the maximum values to try
                    for c0, c1 and c2. (c0 not included if it is fixed)
    initial :       list, Initial guesses for each of the parameters. If c) is
                    fixed, then the initial guess will be used as the fixed 
                    value. [c0, c1, c2]
    test_results :  pd.DF, A dataframe containing the test results from the
                    the recombination tests that were run. 
    
    Returns
    -------------
    [c0, c1, c2] :  list, the estimated coefficients (c0 will be fixed value
                    if it was fixed)
    fit_data :      pd.DF, dataframe containing the relevant fitted values
                    also contains the fitted values for the neher paper

    """
    #Check that we received the correct amount of bounds
    if c0_fixed:
        correctLen = 2
    else: correctLen = 3

    if len(upper_bounds)!= correctLen or len(lower_bounds) != correctLen:
        raise ValueError('Too many or too few bounds given')    

    #Run the fit
    if c0_fixed:
        c0_val = c0_fixed
        c1_c2_list = initial
        print(lower_bounds)
        print(upper_bounds)
        print(c0_val)
        print(c1_c2_list)
        res_lsq = optimize.least_squares(fun = residual, x0 = c1_c2_list, 
                                        bounds = [lower_bounds, upper_bounds], 
                                        kwargs={'dDeltaT' : test_results['window'],
                                                'data': test_results['recomb_frequencies'],
                                                'rec_error' : test_results['Recomb Error'],
                                                'fixed_c0' : c0_val}) 
        c0_estimate = c0_val
        c1_estimate = res_lsq.x[0]
        c2_estimate = res_lsq.x[1]
        print(c0_estimate)
        print(c1_estimate)
        print(c2_estimate)
    else:
        res_lsq = optimize.least_squares(fun = residual, x0 = initial, 
                                        bounds = [lower_bounds, upper_bounds], 
                                        kwargs={'dDeltaT' : test_results['window'],
                                                'data': test_results['recomb_frequencies'],
                                                'rec_error' : test_results['Recomb Error']}) 
        c0_estimate = res_lsq.x[0]
        c1_estimate = res_lsq.x[1]
        c2_estimate = res_lsq.x[2]      
    
    #Make the dataframe of fitted values
    x_vals = list(range(0, max(test_results['window'])))
    fitted_vals = [neher_leitner(x, c0_estimate, c1_estimate, c2_estimate) for x in x_vals]
    fitted_vals_paper = [neher_leitner(x, 0.1, 0.26, .0000439) for x in x_vals]    
    fit_data = pd.DataFrame(list(zip(x_vals, fitted_vals, fitted_vals_paper)), 
                    columns= ['x_vals', 'fitted_vals', 'fitted_vals_paper'])
    
    return [c0_estimate, c1_estimate, c2_estimate], fit_data

###################### Functionality for Bootstrapping ########################
def bootstrap_rho(d_ratio_df, num_boots):
    """ Takes a d_ratio dataframe and resamples from it to bootstrap a rho
    estimate.
    ---------------------------------------------------------------------------
    Params
    ------------
    d_ratio_df: pd.DataFrame, containing the d' ratio for a pair of loci
                also includes d' at the first time point, and information 
                about the two timepoints of sampling
    num_boots:  int, the number of bootstrap resamples to perform

    Returns
    -------------
    lower_fit:  pd.Series, row of the estimate dataframe containing the fit at
                the 2.5 percentile
    upper_fit:  pd.Series, row of the estimate dataframe containing the fit at 
                the 97.5 percentile
    estimate_df:pd.DataFrame, the dataframe containing all of the bootstrapped
                fits in the distribution
    """
    #make a list of all the segregating loci
    all_seg_loc_1 = set(d_ratio_df["Locus_1"].unique())
    all_seg_loc_2 = set(d_ratio_df["Locus_2"].unique())
    all_seg_loc = all_seg_loc_1.union(all_seg_loc_2)
    all_seg_loc = np.array(list(all_seg_loc))
    num_seg = len(all_seg_loc)

    estimates = []

    #sample with replacement num_boots times
    for i in range(num_boots):
        #sample a given number of segregating loci
        seg_loc_sample = np.random.choice(all_seg_loc, size = num_seg, replace = True)
        seg_loc_sample = set(seg_loc_sample)

        #get only the autocorrelation of the chosen loci
        stat_df_sample = d_ratio_df[d_ratio_df["Locus_1"].isin(seg_loc_sample)]
        stat_df_sample = stat_df_sample[stat_df_sample["Locus_2"].isin(seg_loc_sample)]

        coeffs, fit_dat = optimize.curve_fit(neher_leitner, stat_df_sample['Dist_X_Time'], stat_df_sample['d_ratio'], p0 = [0, 0.26, .0000439], maxfev = 10000)
        estimates.append([coeffs[0], coeffs[1], coeffs[2], coeffs[1] * coeffs[2]])
    estimate_df = pd.DataFrame(estimates, columns = ['c0', 'c1', 'c2', 'Estimated_Rho'])
    conf_int = (np.quantile(estimate_df['Estimated_Rho'], 0.025), np.quantile(estimate_df['Estimated_Rho'], 0.975))

    #get the row at the given quantile
    lower_fit = estimate_df.iloc[(estimate_df['Estimated_Rho']-conf_int[0]).abs().argsort()[:2]]
    lower_fit = lower_fit.head(1)
    upper_fit = estimate_df.iloc[(estimate_df['Estimated_Rho']-conf_int[1]).abs().argsort()[:2]]
    upper_fit = upper_fit.head(1)
    return lower_fit, upper_fit, estimate_df