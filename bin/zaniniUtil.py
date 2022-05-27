import pandas as pd
import numpy as np
import sys
import os

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
                
def filter_genotype_df(genotypeDF, segregatingLoci, allele_cutoff,  hap_cutoff):
    """Takes a dataframe of genotypes and filter it based on frequency of the 
    alleles present
    Params
    ------------
    genotypeDF : pd.dataframe, dataframe with two columns indicating loci pairs
                 an aditional column indicates the 2 haplotype at those loci 
                 and the timepoint is also labeled. Data is from one patient 
    segregatingLoci : pd.dataframe, dataframe with the position of each segregating
                      site as well as all the alleles and their frequencies
    allele_cutoff : float, the frequency each allele in a haplotype has to reach
                    for it to be included in the output
    hap_cutoff : float, the frequency each haplotype has to reach for it to be
                 included in the output.
    """
    filtered_genotypes = []
    timepoint_list = segregatingLoci['timepoint'].unique()
    for curr_timepoint in timepoint_list:
        curr_seg = segregatingLoci[segregatingLoci['timepoint'] == curr_timepoint]
        curr_gen = genotypeDF[genotypeDF['timepoint'] == curr_timepoint]

        #put all of the allele frequencies into a dictionary of dictionaries
        #the first key is the locus and the second key is the allele
        freqDict = {}
        for index, row in curr_seg.iterrows():
            locus = row.get('position')
            allele_dict = {}
            allele_dict[row.get('allele_1')] = row.get('freq_1')
            allele_dict[row.get('allele_2')] = row.get('freq_2')
            allele_dict[row.get('allele_3')] = row.get('freq_3')
            allele_dict[row.get('allele_4')] = row.get('freq_4')
            freqDict[locus] = allele_dict


        #now use the dictionary to check the frequency of each allele is above
        #the cutoff
        for index, row in curr_gen.iterrows():
            #get the loci to check
            locus1 = row.get("Locus_1")
            locus2 = row.get("Locus_2")
            haplotype = row.get("2_Haplotype")
            allele1 = haplotype[0]
            allele2 = haplotype[1]

            #get the allele frequencies
            check_1 = (freqDict[locus1])[allele1]
            check_2 = (freqDict[locus2])[allele2]
            if check_1 > allele_cutoff and check_2 > allele_cutoff:
                if check_1 < 0.5:
                    check_1 = 1 - check_1
                if check_2 < 0.5:
                    check_2 = 1-check_2
                row['pA'] = check_1
                row['pB'] = check_2
                filtered_genotypes.append(row)
    filtered_genotypes = pd.DataFrame(filtered_genotypes)

    second_filtered = []
    #Now we want to make sure all our genotypes have frequency > hap_cutoff
    for name, group in filtered_genotypes.groupby(['Locus_1', 'Locus_2', 'timepoint']):
        group_sum = group['Count'].sum()
        group['hap_freq'] = group['Count']/group_sum
        second_filtered.append(group[group['hap_freq'].gt(hap_cutoff)])
    
    filtered_genotypes = pd.concat(second_filtered)
    return filtered_genotypes

def make_viral_load_df(viralLoadDir):
    """ Takes in the path to a folder holding the Zanini viral load data. Reads
    in this data and returns a pandas dataframe containing it.
    ---------------------------------------------------------------------------
    Params
    ------
    viralLoadDir:       str, a string with the path to the directory holding
                        the files with the viral load data
    Returns
    -------
    viralLoadData:      pd.DataFrame, its columns are 'Days from infection',
                        'Viral load [virions/ml]', and 'Participant
    """
    available_files = os.listdir(viralLoadDir)

    #data frame of viral loads for each patient
    viralLoadData = []

    #loop through the files to add to our dataframe
    for curr_file in available_files:
        if curr_file[0] == '.':
            continue
        curr_par = curr_file.split('_')[1]
        curr_par = curr_par.split('.')[0]
        
        #get the current viral loads
        curr_vls = pd.read_csv(viralLoadDir + curr_file, sep = '\t')
        curr_vls['Participant'] = curr_par

        viralLoadData.append(curr_vls)

    viralLoadData = pd.concat(viralLoadData, ignore_index= True)

    return viralLoadData