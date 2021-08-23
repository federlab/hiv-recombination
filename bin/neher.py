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
    """
    groupedGenDF = genotypeDF.groupby(['Locus_1', 'Locus_2'])
    
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

        #loop to -1 since last timepoint has nothing following it
        for curr_time in time_list[:-1]:
            #this is all the haplotypes at this timepoint
            curr_haps_df = group[group['timepoint'] == curr_time]
            #now we need to check if there are three haplotypes that we can use
            #to see recombination
            check_three_haps(curr_haps_df['2_Haplotype'].to_list())
            print(curr_haps_df, file = sys.stderr)
            break
        break

    return

############################ Helper Functions #################################

testCase = 

def check_three_haps(hap_list):
    """ Takes in a list of haplotypes and returns false if there are not three 
    haplotypes that we can use for the haplotype test. Otherwise, returns a 
    list of haplotypes that would satisfy the test
    """
    #make a list to add haplotypes to that would satisfy the test
    return_haps = []

    #position 1
    first_alleles = np.array([x[0] for x in hap_list])
    second_alleles = np.array([x[1] for x in hap_list])
    first_alleles_un = np.unique(first_alleles)
    for al in first_alleles_un:
        #check for alleles at this position that have at least two entries
        #aka are present in two genotypes
        if first_alleles.count(al) < 2:
            continue
        
        #get the alleles its matched with
        matches = [x[1] for x in hap_list if x[0] == al]
        #and alleles we havent observed it with.
        non_matches = []
        for curr_candidate in np.unique(second_alleles):
            if al + curr_candidate not in hap_list:
                non_matches.append(curr_candidate)
        
        #make sure there are some haplotypes we havent observed
        if non_matches == []:
            continue

        #now check the alleles its matched with to see if they have at least
        #two entries
        locus_2_allele1 = []
        locus_1_allele2 = []
        for curr_match in matches:
            if second_alleles.count(curr_match) < 2:
                continue
            else: 
                #put in options for allele at second site
                locus_2_allele1.append(curr_match)
                #store matches corresponding with this allele at the first site
                locus_1_allele2.extend([x[0] for x in hap_list if x[1] == curr_match])
        
        #check if there are options for an allele at the second site
        if locus_2_allele1 == []:
            continue
        #if there are, get possible recombination haplotypes
        for curr_test in locus_2_allele1:




    #position 2
    second_alleles = np.array([x[1] for x in hap_list])
    second_alleles = np.unique(second_alleles)
    print(first_alleles, file = sys.stderr)
    print(second_alleles, file = sys.stderr)
    
    #which position we are checking
    positions = range(2)
    for i in positions:



        



# def enumerateHaplotypes(wideGenDF, locus1, locus2):
#     """Takes a wide dataframe of genotypes and enumerates all the haplotypes
#     across two loci. Returns a Dataframe where each row is a genotype and a
#     count column keeps track of the times each was observed"""
#     #group our dataframe by pairs of loci
#     haplotypesDF = wideGenDF.groupby([locus1,locus2]).size()
#     haplotypesDF.reset_index()
#     haplotypesDF.rename(columns={0:'count'})
#     return haplotypesDF

#first we need to start by looking for three different haplotypes at each locus
# def three_hap_loci(genDF, locus1, locus2):
#     """ Start with a genotype dataframe for a specific patient at a specific
#     time point. 
#     """
#     #convert our dataframe to wide
#     wideGenDF = genDF.pivot(index = 1,columns = 0, values = 2).reset_index()

#     #get the rows where both loci are present
#     currLoc1 = wideGenDF[wideGenDF[locus1].notnull()]
#     bothLoc = currLoc1[currLoc1[locus2].notnull()]

#     #if there are no reads with genotypes at both loci
#     if bothLoc.empty:
#         return False
    
#     #now check that there are three haplotypes
#     #need to make sure only two alleles present at each locus
#     alleles1 = wideGenDf[locus]
#     uniqueAlleles1 = alleles1.unique()
#     alleles2 = wideGenDF[locus]
#     uniqueAlleles2 = alleles2.unique()


#     if len(uniqueAlleles1) == 2 and len(uniqueAlleles2) == 2:
        

#mutation - one time point, one bi-allelic and second mono-allelic. at following time point, both bi-allelic
#A-C, T-A -> AA or TC, recombination evidence, AT - mutation evidence

#get a list of all sites at least bi-allelic at some time point
#current status
#dataframes with allele at each read and time point, 
#list of all genotype pairs

#mark time point on column

#for each patient
    #for each time point from 1 to n-1
        #for each locus1
            #for each locus2>locus1
                #enumerate all genotypes, g1
                #at time point tp+1
                #enumerate all genotypes, g2
                #if list g1 has 3 haplotypes, 
                    #return l1, l2, t1,t2, bool (does g2 have a fourth haplotype)
                #if g1 has 2 or 3 haplotypes,
                    #does g2 have a haplotype that can be created by mutation but not recombination


#work locus by locus, one locus, compare to all other loci that have intersecting reads
#keep track of patient id in each dataframe, also keep track of number of reads supporting each locus specifically
#low frequency variant filtering, possible. 
#also mark down number of reads in each haplotype - 
#edge cases

#AT - AT
#AC - AC
#TC - TC
#     AA
#are some time points having loci covered by vastly different numbers of reads

    
    # for locus1 in lociList:
    #     #get dataframe for our current locus
    #     currLoc1 = wideGenDF[wideGenDF[locus1].notnull()]
    #     for locus2 in lociList:
    #         #order doesn't matter
    #         if locus1 > locus2: