import sys
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
    mutation_df :   
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

        #loop to -1 since last timepoint has nothing following it
        for curr_time in time_list[:-1]:
            #this is all the haplotypes at this timepoint
            curr_haps_df = group[group['timepoint'] == curr_time]

            first_time_haps = set()
            #see if there are any haplotypes we haven't observed yet
            for curr_hap in curr_haps_df['2_Haplotype']:
                if curr_hap not in observed_haps:
                    first_time_haps.add(curr_hap)

            #do we need to check for the 4th haplotype?
            if passed_3_haps is not False:
                #check for the fourth haplotype and mark its distance
                for hap_four in passed_3_haps:
                    if hap_four in first_time_haps:
                        print(name, file = sys.stderr)
                        # recombination_df.append([name])
                #mark the distance between the loci


            #now we need to check if there are three haplotypes that we can use
            #to see recombination on the next pass
            passed_3_haps = check_three_haps(curr_haps_df['2_Haplotype'].to_list())

            #update the haplotype list with everything we observed
            observed_haps.extend(curr_haps_df['2_Haplotype'].to_list())
            break
        break

    return

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
