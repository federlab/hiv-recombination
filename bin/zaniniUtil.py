import pandas as pd
import numpy as np

def find_segregating_diagonal(coCounts_arr):
    """ Takes a numpy array of co-SNP counts from the Zanini et al. data.
    Scans the diagonal of the array to find segregating sites. 
    returns
    --------
    segregatingLoci:    pd.df, dataframe with the position of each segregating
                        site as well as the two most frequent alleles.
    """
    segregatingLoci = []

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
                segregatingLoci.append([i, allele_freqs[0][0], allele_freqs[0][1],
                                    allele_freqs[1][0], allele_freqs[1][1]])

    segregatingLoci = pd.DataFrame(segregatingLoci,
                            columns= ["position", "allele_1", "freq_1",
                                    "allele_2", "freq_2"])
    return segregatingLoci