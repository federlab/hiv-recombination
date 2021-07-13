import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
import os
import r2Analysis as r2
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

fragStarts = {}
fragStarts['p5'] = {'F1': 1, 'F2' : 1469, 'F3' : 3111, 'F4' : 4376, 
                  'F5' : 5424, 'F6' : 7218}

#This file is primarliy for testing our Zanini R^2 pipeline
#I want to print the loci and the cocounts I'm pulling to make sure the
#indexing matches
dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/'
snpPairs = dataDir + "snpPairs/" 
snpCalls = dataDir + "snpCalls/"  
outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/zanini/'

#Get the loci from our snv file
#I chose a file with a decent number of segregating sites
snp_loci = r2.get_Zanini_SNP_Loci(snpCalls + "p5/1057_days.tsv")

#filter our loci so they are only segregating (allele frequency > 1%)
snp_loci = snp_loci[snp_loci["num_seg"] > 1]
snp_loci["sample_t"] = "1057"
snp_loci["participant"] = "p5"

#load our counts array
rel_file = snpPairs + "cocounts_p5_sample_4_F5.npy"
coCounts_arr = np.load(rel_file)

#get the positions spanned by the fragment
fragmentLen = coCounts_arr.shape[-1]
currStart = fragStarts['p5']['F5']
currEnd = currStart + fragmentLen
print("fragment start = " + str(currStart), file = sys.stderr)
print("segregating loci are : ", file = sys.stderr)

#subset the loci to get only those in the fragment
snp_loci = snp_loci.loc[snp_loci.index >= currStart-1]
snp_loci = snp_loci.loc[snp_loci.index <= currEnd]

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

r2List, distList, supportList = r2.calculate_R2_pairCounts(coCounts_arr, segregatingLoci, verbose = True)
r2_results = pd.DataFrame(r2List, columns=['r2'])
r2_results['dist'] = distList
r2_results['support'] = supportList

print(r2_results, file = sys.stderr)

#plot the results for the current participant
sns.set(rc={'figure.figsize':(15,5)})
myplot = sns.scatterplot(x = 'dist', y = 'r2', data = r2_results)
plt.ylim(-0.1,1.1)
plt.xlim(-10,max(r2_results['dist']))
plt.xlabel("Distance Between Loci")
plt.ylabel("R^2 Value")
plt.tight_layout()
plt.savefig(outDir + "p5f5Test")
plt.close()