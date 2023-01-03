# import pysam
# import allel
import sys
import numpy as np
import pandas as pd
import math


################## Functions Written to Work with Slim Dataset ################

def calculate_R2_df(genotype_df, segregating_loci, verbose = False, 
                    statistic = 'r2', saveData = False):
    """ Takes an array of snps at specific loci with counts of the number of 
    times they were observed together. Calculates the R^2 value for each pair.
    ---------------------------------------------------------------------------
    Params
    ---------
    genotypeDF : pd.dataframe, dataframe with two columns indicating loci pairs
                an aditional column indicates the 2 haplotype at those loci 
                and the timepoint is also labeled.
    segregating_loci : a dataframe indicating the segregating loci and the 
                        highest frequency alleles at those loci
    verbose: bool, whether to print out information for debugging
    statistic: str, 'r2', 'D_prime', or 'D' to indicate which statistic to return
    saveData: bool, whether to save our results to a dataframe
    """
    stat_list = []
    distList = []
    supportList = []  
    if saveData:
        resultsDf = []  

    #Add columns for the alleles to our segregating loci dataframe
    genotype_df['Allele_1'] = [x[0] for x in genotype_df['2_Haplotype']]  
    genotype_df['Allele_2'] = [x[1] for x in genotype_df['2_Haplotype']]  

    #Loop through our pairs of loci
    for locus_1 in segregating_loci.index:
        #get the indices for the segregating alleles at this locus
        locus_1_entry = segregating_loci.iloc[[locus_1]]
        i = locus_1_entry["position"].tolist()[0]
        i_allele1 = locus_1_entry['allele_1'].tolist()[0]
        i_allele2 = locus_1_entry['allele_2'].tolist()[0]

        for locus_2 in segregating_loci.index:
            #only consider each pair once
            if locus_1 < locus_2:
                #get the indices for the segregating alleles at this locus
                locus_2_entry = segregating_loci.iloc[[locus_2]]
                j = locus_2_entry["position"].tolist()[0]
                j_allele1 = locus_2_entry['allele_1'].tolist()[0]
                j_allele2 = locus_2_entry['allele_2'].tolist()[0]


                #now we get just the entries corresponding to the haplotypes
                currentCounts = genotype_df[genotype_df['Locus_1'] == i]
                currentCounts = currentCounts[currentCounts['Locus_2'] == j]

                AB_obs = currentCounts[(currentCounts['Allele_1'] == i_allele1) & (currentCounts['Allele_2'] == j_allele1)]
                if len(AB_obs) == 0:
                    AB_obs = 0
                else: AB_obs = AB_obs['Count'].tolist()[0]

                Ab_obs = currentCounts[(currentCounts['Allele_1'] == i_allele1) & (currentCounts['Allele_2'] == j_allele2)]
                if len(Ab_obs) == 0:
                    Ab_obs = 0
                else: Ab_obs = Ab_obs['Count'].tolist()[0]

                aB_obs = currentCounts[(currentCounts['Allele_1'] == i_allele2) & (currentCounts['Allele_2'] == j_allele1)]
                if len(aB_obs) == 0:
                    aB_obs = 0
                else: aB_obs = aB_obs['Count'].tolist()[0]

                ab_obs = currentCounts[(currentCounts['Allele_1'] == i_allele2) & (currentCounts['Allele_2'] == j_allele2)]
                if len(ab_obs) == 0:
                    ab_obs = 0
                else: ab_obs = ab_obs['Count'].tolist()[0]

                allSum = AB_obs + aB_obs + Ab_obs + ab_obs
                #make sure there are enough spanning reads
                if allSum < 10:
                    continue
                #calculate both statistics if we are saving things
                if saveData:
                    r_squared = r2(AB_obs=AB_obs, Ab_obs=Ab_obs,
                                 aB_obs=aB_obs, ab_obs= ab_obs)
                    d_prime = calc_D_prime(AB_obs=AB_obs, Ab_obs=Ab_obs,
                                 aB_obs=aB_obs, ab_obs= ab_obs)
                    d_val = calc_D(AB_obs=AB_obs, Ab_obs=Ab_obs,
                                 aB_obs=aB_obs, ab_obs= ab_obs)
                    p_A = (Ab_obs + AB_obs)/float(allSum)
                    p_B = (AB_obs + aB_obs)/float(allSum)
                    dataList = [i, j, p_A, p_B, AB_obs, Ab_obs, aB_obs, ab_obs, r_squared, d_prime, d_val ]
                    resultsDf.append(dataList)

                if statistic == 'D_prime':
                    curr_stat = calc_D_prime(AB_obs=AB_obs, Ab_obs=Ab_obs,
                                 aB_obs=aB_obs, ab_obs= ab_obs)
                elif statistic == 'D':
                    curr_stat = calc_D(AB_obs=AB_obs, Ab_obs=Ab_obs,
                                 aB_obs=aB_obs, ab_obs= ab_obs)
                else:
                    curr_stat = r2(AB_obs=AB_obs, Ab_obs=Ab_obs,
                                 aB_obs=aB_obs, ab_obs= ab_obs)
                #If the frequencies were too low
                if curr_stat is None:
                    continue

                stat_list.append(curr_stat)
                distance = j - i 
                distList.append(distance)
                support = allSum
                supportList.append(support)
                
                if verbose:
                    print("--------------------")
                    print("Locus 1 is : " + str(i) + " Locus 2 is : " + str(j))
                    print(curr_stat)
                    print("AB")
                    print(AB_obs/allSum)
                    print("aB")
                    print(aB_obs/allSum)
                    print("Ab")
                    print(Ab_obs/allSum)
                    print("ab")
                    print(ab_obs/allSum)

    if saveData:
        resultsDf = pd.DataFrame(resultsDf, columns=[
            'Locus_1', 'Locus_2', 'p_A', 'p_B', 'AB_obs', 'Ab_obs', 
            'aB_obs', 'ab_obs', 'r_squared', 'd_prime', 'd_val'])
        return stat_list, distList, supportList, resultsDf
    return stat_list, distList, supportList

################# Functions Written to Work with Zanini Dataset ###############

def calculate_R2_pairCounts(coCounts_arr, segregating_loci, verbose = False, statistic = 'r2',
                            saveData = False):
    """ Takes an array of snps at specific loci with counts of the number of 
    times they were observed together. Calculates the R^2 value for each pair.
    ---------------------------------------------------------------------------
    Params
    ---------
    coCounts_arr : array of dimension 6 x 6 x L x L where six is the length of
                   the alphabet "ACGT-N" and the last two dimensions refer to
                   the positions of the two alleles
    segregating_loci : a dataframe indicating the segregating loci and the 
                        highest frequency alleles at those loci
    verbose: bool, whether to print out information for debugging
    statistic: str, 'r2' or 'D' to indicate which statistic to return
    saveData: bool, whether to save our results to a dataframe
    """
    stat_list = []
    distList = []
    supportList = []  
    if saveData:
        resultsDf = []  
    #dictionary indicating which alleles correspond to which positions
    #in the co_Counts array
    alleleToPos = {'A' : 0, 'C' : 1, 'G' : 2, 'T' : 3}

    #Loop through our pairs of loci
    for locus_1 in segregating_loci.index:
        #get the indices for the segregating alleles at this locus
        locus_1_entry = segregating_loci.iloc[[locus_1]]
        i = locus_1_entry["position"].tolist()[0]
        i_allele1 = locus_1_entry['allele_1'].tolist()[0]
        i_allele1 = alleleToPos[i_allele1]
        i_allele2 = locus_1_entry['allele_2'].tolist()[0]
        i_allele2 = alleleToPos[i_allele2]
        for locus_2 in segregating_loci.index:
            #only consider each pair once
            if locus_1 < locus_2:
                #get the indices for the segregating alleles at this locus
                locus_2_entry = segregating_loci.iloc[[locus_2]]
                j = locus_2_entry["position"].tolist()[0]
                j_allele1 = locus_2_entry['allele_1'].tolist()[0]
                j_allele1 = alleleToPos[j_allele1]
                j_allele2 = locus_2_entry['allele_2'].tolist()[0]
                j_allele2 = alleleToPos[j_allele2]

                #now we get just the entries corresponding to the haplotypes
                currentCounts = coCounts_arr[:,:,i ,j]
                AB_obs = int(currentCounts[i_allele1, j_allele1])
                Ab_obs = int(currentCounts[i_allele1, j_allele2])
                aB_obs = int(currentCounts[i_allele2, j_allele1])
                ab_obs = int(currentCounts[i_allele2, j_allele2])

                allSum = AB_obs + aB_obs + Ab_obs + ab_obs
                #make sure there are enough spanning reads
                if allSum < 10:
                    continue
                #calculate all three statistics if we are saving things
                if saveData:
                    r_squared = r2(AB_obs=AB_obs, Ab_obs=Ab_obs,
                                 aB_obs=aB_obs, ab_obs= ab_obs)
                    d_prime = calc_D_prime(AB_obs=AB_obs, Ab_obs=Ab_obs,
                                 aB_obs=aB_obs, ab_obs= ab_obs)
                    d_val = calc_D(AB_obs=AB_obs, Ab_obs=Ab_obs,
                                 aB_obs=aB_obs, ab_obs= ab_obs)
                    p_A = (Ab_obs + AB_obs)/float(allSum)
                    p_B = (AB_obs + aB_obs)/float(allSum)
                    dataList = [i, j, p_A, p_B, AB_obs, Ab_obs, aB_obs, ab_obs, r_squared, d_prime, d_val]
                    resultsDf.append(dataList)

                else:
                    if statistic == 'D_prime':
                        curr_stat = calc_D_prime(AB_obs=AB_obs, Ab_obs=Ab_obs,
                                    aB_obs=aB_obs, ab_obs= ab_obs)
                    elif statistic == 'D':
                        curr_stat = calc_D(AB_obs=AB_obs, Ab_obs=Ab_obs,
                                    aB_obs=aB_obs, ab_obs= ab_obs)
                    else:
                        curr_stat = r2(AB_obs=AB_obs, Ab_obs=Ab_obs,
                                    aB_obs=aB_obs, ab_obs= ab_obs)
                    #If the frequencies were too low
                    if curr_stat is None:
                        continue

                stat_list.append(curr_stat)
                distance = j - i 
                distList.append(distance)
                support = allSum
                supportList.append(support)
                
                if verbose:
                    print("--------------------")
                    print("Locus 1 is : " + str(i) + " Locus 2 is : " + str(j))
                    print(curr_stat)
                    print("AB")
                    print(AB_obs/allSum)
                    print("aB")
                    print(aB_obs/allSum)
                    print("Ab")
                    print(Ab_obs/allSum)
                    print("ab")
                    print(ab_obs/allSum)
    if saveData:
        resultsDf = pd.DataFrame(resultsDf, columns=[
            'Locus_1', 'Locus_2', 'p_A', 'p_B', 'AB_obs', 'Ab_obs', 
            'aB_obs', 'ab_obs', 'r_squared', 'd_prime', 'd_val'])
        return stat_list, distList, supportList, resultsDf
    return stat_list, distList, supportList

######################### Calculation Functions ###############################

def r2(AB_obs, Ab_obs, aB_obs, ab_obs):
    """ Takes in the numbers of observations of each haplotype and outputs
    the R^2 statistic
    ---------------------------------------------------------------------------
    returns
    -------
    r_squared: float, The R^2 statistic
    None if frequencies of alleles A or B are greater than 99%
    """
    #sum our haplotypes to get total number of observations
    allSum = AB_obs + aB_obs + Ab_obs + ab_obs

    #calculate the frequencies of our alleles
    p_A = (Ab_obs + AB_obs)/float(allSum)
    p_B = (AB_obs + aB_obs)/float(allSum)

    #make sure we wont get a value error
    if p_A > 0.99 or p_B > 0.99 or p_A < 0.01 or p_B < 0.01: 
        # print("Frequencies were too low")
        return None
    

    topFrac = AB_obs/ (AB_obs + aB_obs + Ab_obs + ab_obs)
    numerator = topFrac - (p_A * p_B)
    numerator = numerator ** 2
    denominator = p_A * p_B * (1-p_A) * (1-p_B)

    r_squared = numerator / denominator
    return r_squared

def calc_D_prime(AB_obs, Ab_obs, aB_obs, ab_obs):
    """ Takes in the numbers of observations of each haplotype and outputs
    the D' statistic
    ---------------------------------------------------------------------------
    returns
    -------
    d_stat: float, The D' statistic
    None if frequencies of alleles A or B are greater than 99%
    """
    #sum our haplotypes to get total number of observations
    allSum = AB_obs + aB_obs + Ab_obs + ab_obs

    #calculate the frequencies of our alleles
    p_A = (Ab_obs + AB_obs)/float(allSum)
    p_B = (AB_obs + aB_obs)/float(allSum)

    #make sure we wont get a value error
    if p_A > 0.99 or p_B > 0.99 or p_A < 0.01 or p_B < 0.01: 
        # print("Frequencies were too low")
        return None
    
    AB_freq = AB_obs/allSum
    Ab_freq = Ab_obs/allSum
    aB_freq = aB_obs/allSum
    ab_freq = ab_obs/allSum
    
    d_stat = AB_freq * ab_freq - aB_freq * Ab_freq
    D_max = None
    if d_stat > 0:
        D_max = min(p_A * (1-p_B), (1-p_A) * p_B)
    elif d_stat == 0:
        return 0
    else:
        D_max = min(p_A * p_B, (1-p_A) * (1-p_B))

    d_stat = d_stat / D_max 
    return abs(d_stat)

def calc_D(AB_obs, Ab_obs, aB_obs, ab_obs):
    """ Takes in the numbers of observations of each haplotype and outputs
    the D statistic. 
    ---------------------------------------------------------------------------
    returns
    -------
    d_stat: float, The D statistic
    None if frequencies of alleles A or B are greater than 99%
    """
    #sum our haplotypes to get total number of observations
    allSum = AB_obs + aB_obs + Ab_obs + ab_obs

    #calculate the frequencies of our alleles
    p_A = (Ab_obs + AB_obs)/float(allSum)
    p_B = (AB_obs + aB_obs)/float(allSum)

    #make sure we wont get a value error
    if p_A > 0.99 or p_B > 0.99 or p_A < 0.01 or p_B < 0.01: 
        # print("Frequencies were too low")
        return None
    
    AB_freq = AB_obs/allSum
    Ab_freq = Ab_obs/allSum
    aB_freq = aB_obs/allSum
    ab_freq = ab_obs/allSum
    
    d_stat = AB_freq * ab_freq - aB_freq * Ab_freq
    return d_stat