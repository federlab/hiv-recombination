import sys
import random
from time import time
import pandas as pd
import numpy as np
from scipy import optimize

#utilize analysis proposed by Neher and Leitner to investigate how viral load
#affects recombination

def run_analysis(genotypeDF, verbose = False, success_filt = 0):
    """Takes a long format dataframe of genotypes at different timepoints. 
    For each pair of loci, 
    ----------------------------------------------------------------------------
    Params
    ------------
    genotypeDF : pd.dataframe, dataframe with two columns indicating loci pairs
                 an aditional column indicates the 2 haplotype at those loci 
                 and the timepoint is also labeled. Data is from one patient 
    verbose : bool, whether to print the tests conducted to a log file. 
    success_filt : a threshold frequency which success haplotypes have to meet
                 for a test to be counted as a pass.

    Returns
    -------------
    recombination_df : pd.dataframe, one column indicates whether a pair of loci
                    triggered a three haplotype test. Two columns give the loci
                    pair and a third column indicates whether they passed the test. 
                    Additionally it has two columns indicating what timepoints these
                    were observed between. The final column indicates the support of
                    the test. 
    mutation_df : pd.dataframe, dataframe similar to recombination df. Indicates
                    when a test has been conducted to see if a new mutation 
                    has been observed. Also indicates the results and support
                    of the test.
    """
    groupedGenDF = genotypeDF.groupby(['Locus_1', 'Locus_2'])
    recombination_df = []
    mutation_df = []
    
    #iterate through pairs of loci
    for name, group in groupedGenDF:
        #get the timepoints for this data
        time_list = group['timepoint'].unique()
        #the labels are strings but they need to be ints to be sorted
        time_list = list(map(int, time_list))
        time_list.sort()

        time_list = list(map(str, time_list))

        #make sure we have multiple timepoints
        if len(time_list) == 1:
            continue

        #whether to check for recombinant haplotypes on next pass
        passed_3_haps = False

        #observed haplotypes
        observed_haps = set()
        #which alleles we've seen at each locus
        observed_alleles_1 = set()
        observed_alleles_2 = set()
        #we dont want any repeat successes
        observed_successes = set()

        #loop to -1 since last timepoint has nothing following it
        for i in range(len(time_list) -1):
            curr_time = time_list[i]

            #this is all the haplotypes at this timepoint
            curr_haps_df = group[group['timepoint'] == curr_time]

            #get the number of reads at the locus at this time point
            supporting_reads = curr_haps_df['Count'].sum()
            p_A = curr_haps_df['pA']
            p_B = curr_haps_df['pB']

            first_time_haps = set()
            #see if there are any haplotypes we haven't observed yet
            for curr_hap in curr_haps_df['2_Haplotype']:
                if curr_hap not in observed_haps:
                    first_time_haps.add(curr_hap)

            #do we need to check for the 4th haplotype (see if there are successes)
            #that have not been previously observed
            if passed_3_haps is not False:
                if verbose:
                    print("***********************\n Possible Successes")
                    print(passed_3_haps)
                    print(curr_haps_df)
                #randomly select a success haplotype
                if len(passed_3_haps) > 1:
                    #get only the successes that haven't been previously observed
                    passed_3_haps = [x for x in passed_3_haps if x not in observed_successes]
                    if len(passed_3_haps) >= 1:
                        passed_3_haps = [random.choice(passed_3_haps)]
                    else: passed_3_haps = False
                     
                #check for the fourth haplotype and mark its distance we dont
                #want any repeat successes either
                if passed_3_haps[0] in first_time_haps and \
                    passed_3_haps[0] not in observed_successes:

                    #compute the frequency of the success haplotype
                    success_freq = curr_haps_df[curr_haps_df['2_Haplotype'] == passed_3_haps[0]]
                    success_freq = success_freq['Count'].tolist()[0]/supporting_reads

                    #check the frequency of the haplotype that triggered the success is high enough
                    if success_freq < success_filt:
                        recombination_df.append([name[0], name[1], False, 
                                        "-", curr_time, time_list[i-1],
                                        supporting_reads, float('inf'), p_A, p_B])
                    else:
                        recombination_df.append([name[0], name[1], True, 
                                        passed_3_haps[0], curr_time, time_list[i-1],
                                        supporting_reads, success_freq, p_A, p_B])
                        observed_successes.add(passed_3_haps[0])
                #save failures so we can calculate frequency
                else:
                    recombination_df.append([name[0], name[1], False, 
                                        "-", curr_time, time_list[i-1],
                                        supporting_reads, float('inf'), p_A, p_B])
                #reset test
                passed_3_haps = False
            
            # # 03/02/2022 Today I moved the mutation test functionality to
            # # The mutation analysis function. I don't think this code bit
            # # -Elena Romero
            # #check for alleles that weren't present before
            # #only do this past the first timepoint though
            # #add all alleles we observe at the first timepoint
            # if i == 0:
            #    for test_hap in first_time_haps:
            #        observed_alleles_1.add(test_hap[0])
            #        observed_alleles_2.add(test_hap[1])
            # else:
            #     mut_found = False
            #     for test_hap in first_time_haps:
            #         allele1 = test_hap[0]
            #         allele2 = test_hap[1]
            #         if allele1 not in observed_alleles_1 \
            #             or allele2 not in observed_alleles_2:
            #             #dont have to worry about duplicates
            #             observed_alleles_1.add(allele1)
            #             observed_alleles_2.add(allele2)
            #             #compute the frequency of the success haplotype
            #             success_freq = curr_haps_df[curr_haps_df['2_Haplotype'] == test_hap]
            #             success_freq = success_freq['Count'].tolist()[0]/supporting_reads

            #             #check the frequency of the haplotype that triggered the success is high enough
            #             if success_freq < success_filt:
            #                 mutation_df.append([name[0], name[1], False, curr_time,
            #                     time_list[i-1], supporting_reads, float('inf'), p_A, p_B]) 
            #             else:
            #                 mutation_df.append([name[0], name[1], True, curr_time,
            #                         time_list[i-1], supporting_reads, success_freq, p_A, p_B])
            #             mut_found = True
            #     #if our test did not find any mutations
            #     if not mut_found:
            #         mutation_df.append([name[0], name[1], False, curr_time,
            #          time_list[i-1], supporting_reads, float('inf'), p_A, p_B])


            #now we need to check if there are three haplotypes that we can use
            #to see recombination on the next pass
            passed_3_haps = check_three_haps(curr_haps_df['2_Haplotype'].to_list())

            #update the haplotype list with everything we observed
            for hap_2 in first_time_haps:
                observed_haps.add(hap_2)

    recombination_df = pd.DataFrame(recombination_df, columns = ["Locus_1", "Locus_2",
                                            "Test_Passed", "Haplotype","Curr_Timepoint", 
                                            "Last_Timepoint", "Supporting_Reads",
                                            "Success_Freq", 'pA', 'pB'])
    return recombination_df

def mutation_analysis(segregatingLoci, fragmentLength, necessary_freq):
    """Takes a long format dataframe of genotypes at different timepoints. 
    Creates a mutation dataframe to calculate the mutation line for the Neher
    and Leitner analysis.
    ----------------------------------------------------------------------------
    Params
    ------------
    segregatingLoci : pd.df, dataframe with the position of each segregating
                 site as well as the frequencies of all alleles at that
                 site.
    fragmentLength : int, the length of the sequenced fragments. This number is
                 incredibly important because it is used to calculate the 
                 number of mutation tests that would be run. (mono-allelic
                 tests are not run explicitly and are counted as failures
                 to save time)
    necessary_freq : float, the frequency necessary for an allele to be counted 
                 as newly observed

    Returns
    -------------
    mutation_df : pd.dataframe, dataframe similar to recombination df. Indicates
                    when a test has been conducted to see if a new mutation 
                    has been observed. Also indicates the results and support
                    of the test. Includes an additional column of test numbers 
                    since the denominator tests are just counted as one entry
                    for multiple tests to save time and space.
    """    
    all_suc_df = []
    no_suc_count = 0
    #get the timepoints to loop throught
    total_time_list = segregatingLoci['timepoint'].unique()

    #a useful set of all four alleles
    nucs = set(['A', 'C', 'G', 'T'])

    #make sure that the times are integers before sorting
    total_time_list = list(map(int, total_time_list))
    total_time_list.sort()
    first_time = total_time_list[0]

    #loop through the loci
    loci_list = segregatingLoci['position']
    for curr_loc in loci_list:
        #get the data for this locus
        curr_seg = segregatingLoci[segregatingLoci['position'] == curr_loc]
        curr_times = curr_seg['timepoint'].unique()
        curr_times = list(map(int, total_time_list))
        curr_times.sort()

        #make a place to store alleles we have observed
        observed_alleles = set()




        #loop through the timepoints
        for i in range(len(curr_times)):
            curr_time = curr_times[i]
            curr_alleles = curr_seg[curr_seg['timepoint'] == curr_time]
            if curr_alleles.empty:
                continue
            present_alleles = \
                get_present_alleles(curr_alleles, necessary_freq)
            new_alleles = present_alleles - observed_alleles

            #check for new alleles if if isn't the first timepoint
            if curr_time != first_time:
                success_allele = nucs - observed_alleles
                if not success_allele:
                    no_suc_count += 1
                    continue

                success_allele = random.choice(list(success_allele))


                if curr_alleles.empty:
                    continue


                #get all the pairwise distances
                distances = list(range(1, curr_loc)) + \
                    list(range(1, fragmentLength - curr_loc))
                
                #make a dataframe of test successes at these distances
                curr_suc_df = pd.DataFrame(distances, columns=['dist'])
                curr_suc_df['Locus'] = curr_loc
                curr_suc_df["Curr_Timepoint"] = curr_time
                curr_suc_df["Last_Timepoint"] = curr_times[i-1]

                #if we are going to have test successes
                if success_allele in new_alleles: 
                    curr_suc_df['Test_Passed'] = True
                    all_suc_df.append(curr_suc_df) 
                    observed_alleles = observed_alleles.union(present_alleles) 
                    break
                else:
                    curr_suc_df['Test_Passed'] = False
                    all_suc_df.append(curr_suc_df)  
                    observed_alleles = observed_alleles.union(present_alleles)
                  
                        
        #add all new alleles to our observations
            observed_alleles = observed_alleles.union(present_alleles)

    #add tests for non-segregating loci to the dataframe
    #loop through the timepoints
    for i in range(len(curr_times)):
        non_seg = set(range(fragmentLength))
        non_seg = non_seg - set(curr_alleles['position'])
        non_seg = list(non_seg)
        for curr_loc in non_seg:
            #get all the pairwise distances
            distances = list(range(1, curr_loc)) + \
                list(range(1, fragmentLength - curr_loc))
            
            #make a dataframe of test successes at these distances
            curr_suc_df = pd.DataFrame(distances, columns=['dist'])
            curr_suc_df['Locus'] = curr_loc
            curr_suc_df['Test_Passed'] = False
            curr_suc_df["Curr_Timepoint"] = curr_time
            curr_suc_df["Last_Timepoint"] = curr_times[i-1]
            all_suc_df.append(curr_suc_df)


    #make a column for distance times time
    all_suc_df = pd.concat(all_suc_df, ignore_index=True) 
    all_suc_df["Dist_x_Time"] = all_suc_df["Curr_Timepoint"] - all_suc_df["Last_Timepoint"]
    all_suc_df["Dist_x_Time"] = all_suc_df["Dist_x_Time"] * all_suc_df['dist']
    print(all_suc_df[all_suc_df['Test_Passed'] == True])


    # #get the denominator
    # all_tests_df = total_mut_tests(fragmentLength, total_time_list)
    # all_tests_df["Dist_x_Time"] = all_tests_df["Curr_Timepoint"] - all_tests_df["Last_Timepoint"]
    # all_tests_df["Dist_x_Time"] = all_tests_df["Dist_x_Time"] * all_tests_df['dist']

    print(no_suc_count)
    #tests are run between each pair of timepoints
    #so once they are returned they need to be multiplied by the num of 
    #timepoints, and the dist column needs to be converted to distance X time
    return all_suc_df #, all_tests_df

##################### Helper Functions for Running Tests ######################
#Test cases for check_three_haps(hap_list)
#return 1 = CG
testCase1 = ['AG', 'AT', 'CT']
#return 2 = False
testCase2 = ['AG', 'CT', 'TC']
#return 3 = [ CG, TT ]
testCase3 = ['AT', 'CT', 'AG', 'TG']
#return 4 = False
testCase4 = ['AT', 'CG', 'CC']
#Already has four haploypes
#return 5 = false
testCase5 = ['AA', 'AT', 'TA', 'TT']
#return 6 = ‘CA’
testCase6 = ['AA', 'AT', 'TA', 'TT', 'CT']

def check_three_haps(hap_list):
    """ Takes in a list of haplotypes and returns false if there are not three 
    haplotypes that we can use for the haplotype test. Otherwise, returns a 
    list of haplotypes that would satisfy the test
    """
    #check that there are at least three haplotypes
    if len(hap_list) < 3:
        return False

    #make a list to add haplotypes to that would satisfy the test
    return_haps = []
    no_overlaps = []

    #get pairs of unique alleles (no overlapping alleles)
    for i in range(len(hap_list)):
        haplotype1 = hap_list[i]
        for j in range(i, len(hap_list)):
            haplotype2 = hap_list[j]
            #check if there are overlapping alleles
            if haplotype1[0] != haplotype2[0] and \
                haplotype1[1] != haplotype2[1]:
                no_overlaps.append((haplotype1, haplotype2))
    
    #now for each of these pairs, check if only combo of the alleles is present
    for curr_pair in no_overlaps:
        combo1 = curr_pair[0][0] + curr_pair[1][1]
        combo2 = curr_pair[1][0] + curr_pair[0][1]
        if combo1 not in hap_list and combo2 in hap_list:
            return_haps.append(combo1)
        if combo2 not in hap_list and combo1 in hap_list:
            return_haps.append(combo2)
    
    if return_haps == []:
        return False
    else: return np.unique(np.array(return_haps))

def estimate_recombination_rate(c0, c1, c2):
    """ Takes in the coefficients from fitting the curve to our data. Then
    returns our estimate of the recombination rate
    ---------------------------------------------------------------------------
    Params
    ------------
    c0 : float, the intercept of our curve fit
    c1 : float, coefficient for second term
    c2 : float, coefficient in the exp

    Returns
    -------------
    rec_rate : the per virus recombination rate
    """
    numerator = c1 * c2 
    denominator = np.log(1/(1-c0))
    denominator = np.log(1/(1 - c0 - c1)) - denominator
    denominator = (1-c0) * denominator
    return numerator/ denominator

def get_present_alleles(allele_freqs, necessary_freq):
    """ Takes in a series in the form of a row from the segregating loci df.
    Returns a set of the alleles that are segregating in that row if they
    are above the necessary frequency
    """
    allele_set = set()
    #loop through the four alleles
    for i in range(1,5):
        curr_freq = allele_freqs['freq_' + str(i)].tolist()[0]
        if curr_freq > necessary_freq:
            allele_set.add(allele_freqs['allele_' + str(i)].tolist()[0])
        #break if the current allele has frequency 0
        #since the frequencies are sorted, the others will all be 0
        else:
            break
    return allele_set

# def total_mut_tests(fragmentLength, time_list):
#     """ Calculates the denominator for the mutation test line by determining
#     the total number of tests that would be run at what distances if we 
#     compared every single pair of loci in a test for one timepoint.
#     ---------------------------------------------------------------------------
#     Params
#     ------------
#     fragmentLength : int, the length of the sequenced fragments. 
#     time_list : list, a list of timepoints (in days or generations) that tests 
#                      will be conducted at. Used to create the distance x time 
#                      column

#     Returns
#     ------------
#     mut_denom_df : pd.df, a pandas dataframe of distances where each distance 
#                 has a given number of tests associated with it.
#     """
#     #A variable for the current number of tests
#     num_tests = 1
#     mut_denoms = []

#     #loop through all possible distances
#     for curr_dist in range(fragmentLength -1, 0, -1):
#         mut_denoms.append([curr_dist, num_tests])
#         num_tests += 1
#     mut_denom_df = pd.DataFrame(mut_denoms, columns= ['dist', 'num_tests'])
#     print(mut_denom_df['num_tests'].sum())
#     time_tests = []
#     #loop through the timepoints and copy the dataframe
#     for i in range(1, len(time_list)):
#         curr_time = time_list[i]
#         curr_df = mut_denom_df.copy()

#         #label the timepoints
#         curr_df['Curr_Timepoint'] = curr_time
#         curr_df['Last_Timepoint'] = time_list[i-1]
#         time_tests.append(curr_df)
    
#     time_tests = pd.concat(time_tests, ignore_index=True)

#     return time_tests


######################### Functionality for Fitting ###########################
def neher_leitner(c0, c1, c2, dDeltaT):
    """The function Neher and Leitner fit to determine recombination rate"""
    return (c0 + (c1 * (1 - np.exp(-c2 * dDeltaT))))

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
    model = neher_leitner(c0, c1, c2, dDeltaT)
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
        c0_val = initial[0]
        c1_c2_list = initial[1:]
        res_lsq = optimize.least_squares(fun = residual, x0 = c1_c2_list, 
                                        bounds = [lower_bounds, upper_bounds], 
                                        kwargs={'dDeltaT' : test_results['window'],
                                                'data': test_results['recomb_frequencies'],
                                                'rec_error' : test_results['Recomb Error'],
                                                'fixed_c0' : c0_val}) 
        c0_estimate = c0_val
        c1_estimate = res_lsq.x[0]
        c2_estimate = res_lsq.x[1]
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
    fitted_vals = [neher_leitner(c0_estimate, c1_estimate, c2_estimate, x) for x in x_vals]
    fitted_vals_paper = [neher_leitner(0.1, 0.26, .0000439, x) for x in x_vals]    
    fit_data = pd.DataFrame(list(zip(x_vals, fitted_vals, fitted_vals_paper)), 
                    columns= ['x_vals', 'fitted_vals', 'fitted_vals_paper'])
    
    return [c0_estimate, c1_estimate, c2_estimate], fit_data
