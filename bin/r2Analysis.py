import pysam
import allel
import sys
import numpy as np
import pandas as pd


def get_VCF_Loci(vcfFile):
    """ Open a vcf file and save the variants in it to a dataframe 
    ---------------------------------------------------------------------------
    Params
    -------
    vcfFile     str, path to the vcf file 
    Returns
    -------
    vcfDat   pd.dataframe, containing each variant for the current dataset
    """
    #load our data into a dataframe
    vcfDat = allel.read_vcf(vcfFile)
    #check for extra alternate alleles
    if vcfDat['variants/ALT'].shape[1] > 3:
        print('More than 3 alternate alleles', file = sys.stderr)

    vcfDat = pd.DataFrame({'POS': vcfDat['variants/POS'], 
            'REF': vcfDat['variants/REF'] ,
             'ALT1': (vcfDat['variants/ALT'])[:,0], 
             'ALT2': (vcfDat['variants/ALT'])[:,1], 
             'ALT3': (vcfDat['variants/ALT'])[:,2]})

    return vcfDat

def get_sam_reads(vcfDF, bamfile):
    """ Loop through a vcf file and get all possible pairs of loci that are
    segregating.
    ---------------------------------------------------------------------------
    Params
    ------
    vcfDF      pd.dataframe, containing each variant for the current dataset
    bamfile    str, the path to the file with the aligned reads
    """
    #get our samfile using pysam
    alignedReads = pysam.AlignmentFile(bamfile, 'rb')

    #make our dataframe of reads
    readGenotypes = []

    #get all of the loci that are needed and make a set
    segregatingLoci = set(vcfDF['POS'])

    #get first and last locus
    firstLocus = vcfDF.iloc[[0]]
    lastLocus = vcfDF.iloc[[-1]]

    #Iterate through the reads
    for pileupcolumn in alignedReads.pileup(start = firstLocus['POS'] , stop = lastLocus['POS'] ):
        #look at all of the loci and record the read's genotype at each locus
        #if the locus is one of the segregating ones
        if pileupcolumn.pos in segregatingLoci:
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    #make a wide row in our dataframe
                    readGenotypes.append([pileupcolumn.pos,
                                pileupread.alignment.query_name,
                                pileupread.alignment.query_sequence[pileupread.query_position] ])
    
    #now make our dataframe
    readGenotypes = pd.DataFrame(readGenotypes)
    #next we need to filter out entries where all of the genotypes are the same
    dfList = []
    #for each locus check how many alleles there are
    for currLoc in (readGenotypes[0]).unique():
        allOccurrences = readGenotypes[readGenotypes[0] == currLoc]
        uniqueAlleles = (allOccurrences[2]).unique()
        #if there is more than one allele at the locus
        if len(uniqueAlleles) > 1:
            dfList.append(allOccurrences)
    filteredGenotypes = pd.concat(dfList)

    return filteredGenotypes
    




def 
    # read  loc1    loc2    loc3
    # 1     #       g       t
    # 2      a       g       #

    #fetch our reads, I think a dataframe ordered by position would be best
    #we need to loop through the reads and for each, record the genotype
    #at all of the segregating positions. We could have a column for each and then
    #do like a pandas loc to find columns where the entry for that locus is not like
    #some prespecified null value



    # #we need to get each pair of loci. To do this we can iterate through the
    # #rows with nested for loops
    # for index1, locus1 in vcfDF.iterrows():
    #     #just get all the reads where locus one is present
    #     for index2, locus2 in vcfDF.iterrows():
    #         #if the two loci are not the same, also we don't care about the
    #         #order of pairs so we can just use greater than to not get 
    #         #the two pairs
    #         if index1 > index2:
                
    #             locusDist = locus1['POS'] - locus2['POS']
    #             #if its possible for the two loci to be on the same read
    #             if locusDist < 201:
    #                 print(locusDist, file = sys.stderr)

                    #locate rows where loc1 not # and loc2 not #

                # for read in alignedReads.fetch(start = locus1['POS'], stop = locus2['POS']):
                #     print(read.is_paired, file = sys.stderr)
                #     print(read.cigartuples, file = sys.stderr)
                #     #now we can use the two loci to get all of the reads between them
                #     #we can do this using pysam
                #     print(locus1, file = sys.stderr)
                #     print(locus2, file = sys.stderr)
                #then we need to check for spanning reads


    return





def calculateR2():
    """
    Given two segregating loci, calculate their R^2 value
    https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0048588
    """
    #We need to know both alleles at each locus
    #To get this we can use our vcf file.


    
    #I think I need to ask more about which files to use (should I use samples
    # from multiple patient dates, or should there be a plot for each date)

    #We need the allele frequencies from the intersecting reads

    