import sys
import random
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
    """
    groupedGenDF = genotypeDF.groupby(['Locus_1', 'Locus_2'])
    recombination_df = []
    
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
            mpApB = supporting_reads * p_A * p_B

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
                    else: 
                        passed_3_haps = False
                        continue

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
            curr_time = str(curr_times[i])
            curr_alleles = curr_seg[curr_seg['timepoint'] == curr_time]
            if curr_alleles.empty:
                continue
            present_alleles = \
                get_present_alleles(curr_alleles, necessary_freq)
            new_alleles = present_alleles - observed_alleles

            #check for new alleles if if isn't the first timepoint
            if int(curr_time) != first_time:
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
                curr_suc_df["Curr_Timepoint"] = int(curr_time)
                curr_suc_df["Last_Timepoint"] = int(curr_times[i-1])

                #if we are going to have test successes
                if success_allele in new_alleles: 
                    curr_suc_df['Test_Passed'] = True
                    all_suc_df.append(curr_suc_df) 
                    observed_alleles = observed_alleles.union(present_alleles) 
                    # break
                else:
                    curr_suc_df['Test_Passed'] = False
                    all_suc_df.append(curr_suc_df)  
                    observed_alleles = observed_alleles.union(present_alleles)
                  
                        
        #add all new alleles to our observations
            observed_alleles = observed_alleles.union(present_alleles)

    #add tests for non-segregating loci to the dataframe
    #loop through the timepoints
    for i in range(1, len(total_time_list)):
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
            curr_suc_df["Curr_Timepoint"] = int(curr_times[i])
            curr_suc_df["Last_Timepoint"] = int(curr_times[i-1])

            all_suc_df.append(curr_suc_df)


    #make a column for distance times time
    all_suc_df = pd.concat(all_suc_df, ignore_index=True) 

    return all_suc_df 

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
