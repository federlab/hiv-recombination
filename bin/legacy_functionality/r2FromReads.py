# import pysam
# import allel
import sys
import numpy as np
import pandas as pd
import math

# def get_VCF_Loci(vcfFile):
#     """ Open a vcf file and save the variants in it to a dataframe 
#     ---------------------------------------------------------------------------
#     Params
#     -------
#     vcfFile     str, path to the vcf file 
#     Returns
#     -------
#     vcfDat   pd.dataframe, containing each variant for the current dataset
#     """
#     #load our data into a dataframe
#     vcfDat = allel.read_vcf(vcfFile)
#     #make sure our dictionary isn't empty
#     if not vcfDat:
#         print('No variants detected in ' + vcfFile, file = sys.stderr)
#         return None

#     vcfDat = pd.DataFrame({'POS': vcfDat['variants/POS'], 
#             'REF': vcfDat['variants/REF'] ,
#              'ALT1': (vcfDat['variants/ALT'])[:,0], 
#              'ALT2': (vcfDat['variants/ALT'])[:,1], 
#              'ALT3': (vcfDat['variants/ALT'])[:,2]})

#     return vcfDat

# def get_sam_reads(vcfDF, bamfile):
#     """ Loop through a vcf file and use it to generate a dataframe of 
#     reads with their genotype at each segregating loci.
#     ---------------------------------------------------------------------------
#     Params
#     ------
#     vcfDF      pd.dataframe, containing each variant for the current dataset
#     bamfile    str, the path to the file with the aligned reads
#     """
#     #get our samfile using pysam
#     alignedReads = pysam.AlignmentFile(bamfile, 'rb')

#     #make our dataframe of reads
#     readGenotypes = []

#     #get all of the loci that are needed and make a set
#     segregatingLoci = set(vcfDF['POS'])

#     #get first and last locus
#     firstLocus = vcfDF.iloc[[0]]
#     lastLocus = vcfDF.iloc[[-1]]

#     #Iterate through the reads
#     for pileupcolumn in alignedReads.pileup(start = firstLocus['POS'] , stop = lastLocus['POS'] ):
#         #look at all of the loci and record the read's genotype at each locus
#         #if the locus is one of the segregating ones
#         if pileupcolumn.pos in segregatingLoci:
#             for pileupread in pileupcolumn.pileups:
#                 if not pileupread.is_del and not pileupread.is_refskip:
#                     print(pileupread.alignment.query_name, file = sys.stderr)
#                     #make a wide row in our dataframe
#                     readGenotypes.append([pileupcolumn.pos,
#                                 pileupread.alignment.query_name,
#                                 pileupread.alignment.query_sequence[pileupread.query_position] ])
    
#     #now make our dataframe
#     readGenotypes = pd.DataFrame(readGenotypes)
#     #next we need to filter out entries where all of the genotypes are the same
#     dfList = []
#     #for each locus check how many alleles there are
#     for currLoc in (readGenotypes[0]).unique():
#         allOccurrences = readGenotypes[readGenotypes[0] == currLoc]
#         uniqueAlleles = (allOccurrences[2]).unique()
#         #if there is more than one allele at the locus
#         if len(uniqueAlleles) > 1:
#             dfList.append(allOccurrences)
#     filteredGenotypes = pd.concat(dfList)

#     return filteredGenotypes


# def r2_all_loci(genDF):
#     """
#     Loop through dataframe of reads and for loci that are present on same read,
#     Calculate R^2 and distance
#     ---------------------------------------------------------------------------
#     Params
#     ------
#     genDF:  pd.DataFrame, each row is a read's genotype at a specific locus
#                             [position, read name, base]
#     Returns:
#     -------
#     r2List: list, R^2 values for each locus
#     distList: list, distances between each corresponding pair of loci
#     supportList: list, number of reads supporting each corresponding r^2 value
#     """
#     #make a list to store the values in 
#     r2List = []
#     distList = []
#     supportList = []

#     #our loci
#     lociList = genDF[0].unique()

#     #convert our dataframe to wide
#     wideGenDF = genDF.pivot(index = 1,columns = 0, values = 2).reset_index()


#     #we need every pair of loci
#     for locus1 in lociList:
#         #get dataframe for our current locus
#         currLoc1 = wideGenDF[wideGenDF[locus1].notnull()]
#         for locus2 in lociList:
#             #order doesn't matter
#             if locus1 > locus2:
#                 currDist = locus1 - locus2
#                 #get all of the reads with a genotype for both loci
#                 bothLoc = currLoc1[currLoc1[locus2].notnull()]
#                 if not bothLoc.empty:
#                     #calculate R^2 for these two loci
#                     currR2 = calculateR2(bothLoc, locus1, locus2)
#                     support = len(bothLoc.index)
#                     #if both the loci are polymorphic for some reads
#                     if currR2 is not None:
#                         r2List.append(currR2)
#                         distList.append(locus1-locus2)
#                         supportList.append(support)
                        
                        
#     return r2List, distList, supportList



# def calculateR2(readDf, locus1, locus2, verbose = False):
#     """
#     Given two segregating loci, calculate their R^2 value using the direct
#     inference method
#     https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0048588
#     Params
#     ------
#     readDF:  pd.DataFrame, each row is a read's genotype at a specific locus
#                         [position, read name, base]. Dataframe should contain
#                         all of the reads that have genotypes at both loci
#     locus1:  int, position of first locus
#     locus2:  int, position of second locus
#     verbose: bool, whether to print debugging info
#     Returns
#     -------
#     None if both loci are not segregating on the overlapping reads
#     r2: float, the calculated R^2 value.
#     """
#     #check how many alleles at the locus
#     numLoc1 = readDf[locus1].unique()
#     numLoc2 = readDf[locus2].unique()

#     #We need to know the frequency of each allele
#     freqList1 = calculate_allele_freq(readDf, locus1)
#     freqList2 = calculate_allele_freq(readDf, locus2)

#     #if both loci are not polymorphic on the overlapping reads, return None
#     if len(freqList1) == 1 or len(freqList2) == 1:
#         return None
    
#     #sort our lists by most frequent
#     freqList1 = sorted(freqList1, key = lambda x: x[1], reverse=True)
#     freqList2 = sorted(freqList2, key = lambda x: x[1], reverse=True)

#     #make sure both our loci have frequency higher than 1%
#     #we just need to check the second entries since the list is sorted
#     if freqList1[1][1] < 0.01 or freqList2[1][1] < 0.01:
#         return None

#     #throw out all alleles that aren't the two most frequent 
#     #and start labeling each genotype
#     freqList1 = freqList1[:2]
#     loc1A = readDf[readDf[locus1] == freqList1[0][0]]
#     loc1a = readDf[readDf[locus1] == freqList1[1][0]]
#     freqList2 = freqList2[:2]

#     #we need to count the observations for each haplotype
#     AB = loc1A[loc1A[locus2] == freqList2[0][0]]
#     aB = loc1a[loc1a[locus2] == freqList2[0][0]]
#     Ab = loc1A[loc1A[locus2] == freqList2[1][0]]
#     ab = loc1a[loc1a[locus2] == freqList2[1][0]]

#     AB_obs = len(AB.index)
#     aB_obs = len(aB.index)
#     Ab_obs = len(Ab.index)
#     ab_obs = len(ab.index)

#     allSum = AB_obs + aB_obs + Ab_obs + ab_obs

#     p_A = (Ab_obs + AB_obs)/float(allSum)
#     p_B = (AB_obs + aB_obs)/float(allSum)

#     #make sure the frequencies aren't to close to 1
#     if p_A == 1 or p_B == 1: return None

#     #now do the calculation
#     topFrac = AB_obs/ (AB_obs + aB_obs + Ab_obs + ab_obs)
#     # print(topFrac, file = sys.stderr)
#     numerator = topFrac - (p_A * p_B)
#     numerator = numerator ** 2
#     # print(numerator, file = sys.stderr)
#     denominator = p_A * p_B * (1-p_A) * (1-p_B)
#     # print(denominator, file = sys.stderr)
    
#     r2 = numerator / denominator

#     if verbose:
#         if r2 < 0.1  or r2 > 0.9:
#             print("--------------------")
#             print("Locus 1 is : " + str(locus1) + " Locus 2 is : " + str(locus2))
#             print(r2)
#             print("AB")
#             print(AB_obs/allSum)
#             print("aB")
#             print(aB_obs/allSum)
#             print("Ab")
#             print(Ab_obs/allSum)
#             print("ab")
#             print(ab_obs/allSum)

#     return r2
    
    
#     #I think I need to ask more about which files to use (should I use samples
#     # from multiple patient dates, or should there be a plot for each date)

#     #We need the allele frequencies from the intersecting reads

# def calculate_allele_freq(readDf, locus): 
#     """Calculate the frequency of an allele at a locus
#     ---------------------------------------------------------------------------
#     Params
#     ------
#     readDF:  pd.DataFrame, each row is a read's genotype at a specific locus
#                         [position, read name, base]. Dataframe should contain
#                         all of the reads that have genotypes at both loci
#     locus:  int, position of locus
    
#     Returns
#     -------
#     freqList:   list of tuples, the first element in each tuple is the allele
#                         and the second element is its frequency. 
#     """
#     #initialize our list to return
#     freqList = []

#     #get how many unique alleles there are
#     currLoc = readDf[locus]
#     uniqueAlleles = currLoc.unique()
#     denominator = len(readDf.index)

#     for allele in uniqueAlleles:
#         readsWithGenotype = readDf[readDf[locus] == allele]
#         count = len(readsWithGenotype.index)
#         frequency = count/float(denominator)
#         freqList.append((allele, frequency))

#     return freqList