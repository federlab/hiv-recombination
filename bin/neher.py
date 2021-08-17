import os
import sys
import numpy as np
import r2Analysis as r2
#utilize analysis proposed by Neher and Leitner to investigate how viral load
#affects recombination

#I need to get started here by putting the data into a dataframe of genotypes.

def plotPoint(genotypeDF, segregatintLoci):
    """Takes a long formats dataframe of genotypes at different timepoints. 
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
        pass
    # #make a dataframe that will be in long format. 
    # resultsDF = []

    # We want to loop through pairs of loci so I think we can just loop through the
    # index of the genotype dataframe
    
    #I think we want to loop through all possible pairs of loci.


    #             currDist = locus1 - locus2
    #             #get all of the reads with a genotype for both loci
    #             bothLoc = currLoc1[currLoc1[locus2].notnull()]
    #             if not bothLoc.empty:
    #                 #now we need to loop through all our timepoints
    #                 for curr_timepoint in timepointList:
    #                     print(bothLoc, file = sys.stderr)
    #                     #get reads at this timepoint
    #                     time_t = bothLoc[bothLoc['timepoint' == curr_timepoint]]
    #                     if time_t.notNull():
    #                         print(time_t, file = sys.stderr)

                        # #enumerate all the haplotypes at this timepoint
                        # print(enumerateHaplotypes(
                        #     bothLoc[bothLoc['timepoint' == curr_timepoint]],
                        #     locus1, locus2), file = sys.stderr)


############################ Helper Functions #################################

def enumerateHaplotypes(wideGenDF, locus1, locus2):
    """Takes a wide dataframe of genotypes and enumerates all the haplotypes
    across two loci. Returns a Dataframe where each row is a genotype and a
    count column keeps track of the times each was observed"""
    #group our dataframe by pairs of loci
    haplotypesDF = wideGenDF.groupby([locus1,locus2]).size()
    haplotypesDF.reset_index()
    haplotypesDF.rename(columns={0:'count'})
    return haplotypesDF

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