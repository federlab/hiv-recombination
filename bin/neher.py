import sys
import pandas as pd
import numpy as np

#utilize analysis proposed by Neher and Leitner to investigate how viral load
#affects recombination

#I need to get started here by putting the data into a dataframe of genotypes.

def run_analysis(genotypeDF, segregatingLoci):
    """Takes a long format dataframe of genotypes at different timepoints. 
    For each pair of loci, 
    Params
    ------------
    genotypeDF : pd.dataframe, dataframe with two columns indicating loci pairs
                 an aditional column indicates the 2 haplotype at those loci 
                 and the timepoint is also labeled. Data is from one patient 
    segregatingLoci : pd.dataframe, dataframe with the position of each segregating
                      site as well as all the alleles and their frequencies
    Returns
    -------------
    recombination_df : pd.dataframe, one column indicates whether a pair of loci
                        triggered a three haplotype test. Two columns give the loci
                        pair and a third column indicates whether they passed the test. 
                        Additionally it has two columns indicating what timepoints these
                        were observed between.
    mutation_df :   
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

        #loop to -1 since last timepoint has nothing following it
        for i in range(len(time_list) -1):
            curr_time = time_list[i]

            #this is all the haplotypes at this timepoint
            curr_haps_df = group[group['timepoint'] == curr_time]

            first_time_haps = set()
            #see if there are any haplotypes we haven't observed yet
            for curr_hap in curr_haps_df['2_Haplotype']:
                if curr_hap not in observed_haps:
                    first_time_haps.add(curr_hap)

            #do we need to check for the 4th haplotype?
            if passed_3_haps is not False:
                #variable to check if successful
                success = False
                #check for the fourth haplotype and mark its distance
                for hap_four in passed_3_haps:
                    if hap_four in first_time_haps:
                        recombination_df.append([name[0], name[1], True, 
                                        hap_four, curr_time, time_list[i-1]])
                        success = True
                #save failures so we can caluclate frequency
                if not success:
                    recombination_df.append([name[0], name[1], False, 
                                        "-", curr_time, time_list[i-1]])
                #reset test
                passed_3_haps = False
            
            #check for alleles that weren't present before
            #only do this past the first timepoint though
            #add all alleles we observe at the first timepoint
            if i == 0:
               for test_hap in first_time_haps:
                   observed_alleles_1.add(test_hap[0])
                   observed_alleles_2.add(test_hap[1])
            else:
                mut_found = False
                for test_hap in first_time_haps:
                    allele1 = test_hap[0]
                    allele2 = test_hap[1]
                    if allele1 not in observed_alleles_1 \
                        or allele2 not in observed_alleles_2:
                        #dont have to worry about duplicates
                        observed_alleles_1.add(allele1)
                        observed_alleles_2.add(allele2)
                        mutation_df.append([name[0], name[1], True, curr_time, time_list[i-1]])
                        mutFount = True
                #if our test did not find any mutations
                if not mut_found:
                    mutation_df.append([name[0], name[1], False, curr_time, time_list[i-1]])


            #now we need to check if there are three haplotypes that we can use
            #to see recombination on the next pass
            passed_3_haps = check_three_haps(curr_haps_df['2_Haplotype'].to_list())

            #update the haplotype list with everything we observed
            for hap_2 in curr_haps_df['2_Haplotype'].to_list():
                observed_haps.add(hap_2)
    recombination_df = pd.DataFrame(recombination_df, columns = ["Locus_1", "Locus_2",
                                            "Test_Passed", "Haplotype","Curr_Timepoint", 
                                            "Last_Timepoint"])
    mutation_df = pd.DataFrame(mutation_df, columns = ["Locus_1", "Locus_2",
                                        "Test_Passed", "Curr_Timepoint", 
                                        "Last_Timepoint"])
    return recombination_df, mutation_df

############################ Helper Functions #################################
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
