import pandas as pd
import numpy as np
import sys

def find_segregating_diagonal(coCounts_arr, all_seg = False):
    """ Takes a numpy array of co-SNP counts from the Zanini et al. data.
    Scans the diagonal of the array to find segregating sites. 
    Params
    ------
    coCounts_arr:   np.array, array from Zanini et al. of snp counts
    all_seg:        bool, if true, all information for all segregating 
                    alleles will be included. This is needed for neher
                    analysis (neher.py).

    returns
    --------
    segregatingLoci:    pd.df, dataframe with the position of each segregating
                        site as well as the two most frequent alleles.
    """
    segregatingLoci = []
    fragmentLen = coCounts_arr.shape[-1]

    #scan our array and get the observed frequency for each row
    for i in range(fragmentLen):
        #I think this is the diagonal which is just counts with itself
        currObs = coCounts_arr[:,:,i,i]
        #!!!!!!! we can use the diagonal to find segregating loci
        #the first item is the a count, the second is the c count and so on
        acgt_counts = np.array(
            [currObs[0][0], currObs[1][1], currObs[2][2], currObs[3][3]])
        
        #If there are multiple alleles check for frequency >1%
        if np.count_nonzero(acgt_counts) > 1:
            total = np.sum(acgt_counts)
            frequencies = acgt_counts/total

            #check if a value is greater than 1%
            checkPercent = lambda x : 1 if (x > 0.01) else 0
            numSegregating = sum(list(map(checkPercent, frequencies)))
            if numSegregating > 1:
                #find the alleles with the highest frequency
                allele_freqs = [("A", frequencies[0]), ("C", frequencies[1]),
                                ("G", frequencies[2]), ("T", frequencies[3])]
                allele_freqs = sorted(allele_freqs, key = lambda x: x[1], reverse=True)
                #append the info for these alleles
                freq_list = [i, allele_freqs[0][0], allele_freqs[0][1],
                                    allele_freqs[1][0], allele_freqs[1][1]]
                #if we want all of them, amend the list
                #this will give us some alleles with zero frequency in the dataframe
                if all_seg:
                    freq_list = [i, 
                                allele_freqs[0][0], allele_freqs[0][1],
                                allele_freqs[1][0], allele_freqs[1][1],
                                allele_freqs[2][0], allele_freqs[2][1],
                                allele_freqs[3][0], allele_freqs[3][1]]
                           
                segregatingLoci.append(freq_list)

    if all_seg: 
        segregatingLoci = pd.DataFrame(segregatingLoci,
                            columns= ["position", "allele_1", "freq_1",
                                    "allele_2", "freq_2", "allele_3", "freq_3",
                                    "allele_4", "freq_4"])
    else: 
        segregatingLoci = pd.DataFrame(segregatingLoci,
                            columns= ["position", "allele_1", "freq_1",
                                    "allele_2", "freq_2"])
    
    return segregatingLoci

def make_genotype_df(segregatingLoci, coCounts_arr):
    """ Takes a list of segregating loci and the of co-SNP counts from the 
    Zanini et al. data. Then makes a dataframe where each row represents 
    a pair of loci and there is an entry containing a genotype observed for 
    that loci pair.
    ---------------------------------------------------------------------------
    Params
    ------
    segregatingLoci:    pd.df, dataframe with the position of each segregating
                        site as well as the two most frequent alleles.
    coCounts_arr:       np.array, array from Zanini et al. of snp counts
    Returns
    -------
    genotype_df:        pd.df, dataframe where each row corresponds to a 
                        2 haplotype observed at a pair of segregating loci.
    """
    #this dictionary matches positions in arrays to alleles
    allele_dict = { 0 : 'A',
                    1 : 'C',
                    2 : 'G',
                    3 : 'T',
                    }
    
    genotype_df = []

    #loop through pairs of segregating loci to get their genotype info
    for locus_1 in segregatingLoci.index:
        locus_1_entry = segregatingLoci.iloc[[locus_1]]
        i = locus_1_entry["position"].tolist()[0]
        for locus_2 in segregatingLoci.index:
            #only consider each pair once
            if locus_1 < locus_2:
                locus_2_entry = segregatingLoci.iloc[[locus_2]]
                j = locus_2_entry["position"].tolist()[0]
                #now pull the array entry for this pair to get the genotypes
                currentCounts = coCounts_arr[:,:,i ,j]
                #from here I think we should loop through the genotype array
                for pos_1 in range(0,4):
                    for pos_2 in range(0,4):
                        gen_count = currentCounts[pos_1,pos_2]
                        #we want to put all the counts that arent zero in our array
                        if gen_count > 0:
                            #but first check that neither of the frequencies are below
                            #the 3% error rate
                            alleles = allele_dict[pos_1] + allele_dict[pos_2]
                            genotype_df.append([i, j, alleles, gen_count])
    
    
    genotype_df = pd.DataFrame(genotype_df, columns= ['Locus_1', 'Locus_2', '2_Haplotype', 'Count'])
    return genotype_df
                
