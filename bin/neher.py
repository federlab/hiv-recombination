import numpy as np
import r2Analysis as r2
#utilize analysis proposed by Neher and Leitner to investigate how viral load
#affects recombination

#Start making background slides
#concepts to touch on 
# why recombination is important in viruses in general
# can increase in frequency and be transmissible, has played an important role in 
# shaping HIV diversity globally
# recombination between two very different variants vs. subtler recombination within patients
# can contribute to immune escape, production of multi-drug resistant variants
# how do we measure recombination rate - past methods
# set up problem : why recombination in hiv might depend on viral load
#
# 
#Speak with Erick
#schedule for some time next week , plots stratified by patient, talk through 
#rough version of slides as far as background
#third check in with erick, what we've gotten done and what we would push to get done
#
#
#Plan to talk about Strauli preprint - maybe early next week. what the data is coming from. 
#Tuesday at 1pm PST, ping erick to chat mid to late next week
#Meet with erick 

#first we need to start by looking for three different haplotypes at each locus
def three_hap_loci(genDF, locus1, locus2):
    """ Start with a genotype dataframe for a specific patient at a specific
    time point. 
    """
    #convert our dataframe to wide
    wideGenDF = genDF.pivot(index = 1,columns = 0, values = 2).reset_index()

    #get the rows where both loci are present
    currLoc1 = wideGenDF[wideGenDF[locus1].notnull()]
    bothLoc = currLoc1[currLoc1[locus2].notnull()]

    #if there are no reads with genotypes at both loci
    if bothLoc.empty:
        return False
    
    #now check that there are three haplotypes
    #need to make sure only two alleles present at each locus
    alleles1 = wideGenDf[locus]
    uniqueAlleles1 = alleles1.unique()
    alleles2 = wideGenDF[locus]
    uniqueAlleles2 = alleles2.unique()


    if len(uniqueAlleles1) == 2 and len(uniqueAlleles2) == 2:
        

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
                    #return l1, l2, t1,t2, bool (does g2 have a haplotype)
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
