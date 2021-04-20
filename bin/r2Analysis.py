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

    return vcfDAt

def get_sam_reads(vcfDF, samFile):
    """ Loop through a vcf file and get all possible pairs of loci that are
    segregating.
    ---------------------------------------------------------------------------
    Params
    ------
    vcfDF      pd.dataframe, containing each variant for the current dataset
    """
    #we need to get each pair of loci. To do this we can iterate through the
    #rows with nested for loops
    for index1, locus1 in vcfDF.iterrows():
        for index2, locus2 in vcfDF.iterrows():
            #now we can use the two loci to get all of the reads between them
            #we can do this using pysam
            print(locus1, file = sys.stderr)
            print(locus2, file = sys.stderr)

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

    