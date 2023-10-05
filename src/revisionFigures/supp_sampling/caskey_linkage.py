import os
import sys
sys.path.append('/Volumes/feder-vol1/home/evromero/2023_hiv-bnabs/bin')
import numpy as np
import pandas as pd
import r2Analysis as r2
from Bio import SeqIO

SEG_CUTOFF = 0.1

#I am going to read in all the data for a single individual and calculate
#the number of segregating sites. This will help me determine if there is
#enough data to proceed with estimation.
###############################################################################
#Read in the data
dataDir = '/Volumes/feder-vol1/home/evromero/2023_hiv-bnabs/data/caskey2017/'
outDir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/caskey/'

#1HB1/1HB1_codon_aligned.fasta

#A dictionary to map timepoints to days
timepoint_dict = {'preinf' : 0,
                  'D0' : 0,
                  'W1' : 7 * 0.5,
                  'W2' : 14 * 0.5,
                  'W4' : 28 * 0.5,
                  'W8' : 56 * 0.5,
                  'W12' : 84 * 0.5,
                  'W16' : 112 * 0.5,
                  'W20' : 140 * 0.5,
                  'W24' : 168 * 0.5}

all_seg_counts = []

#Loop through the directories
for curr_dir in os.listdir(dataDir):
    if curr_dir[0] == '.':
        continue
    if not os.path.isdir(dataDir + curr_dir):
        continue
    inFile = dataDir + curr_dir + '/' + curr_dir + '_codon_aligned.fasta'

    #Read in the fasta sequences
    fasta_sequences = SeqIO.parse(open(inFile),'fasta')

    #Sort the fasta sequences into a dictionary by time point
    fasta_dict = {}
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        #Don't include the reference sequence
        if name == 'B.FR.1983.HXB2-LAI-IIIB-BRU.K03455' or name == 'HXB2':
            continue
        
        #Get the time point
        seq_info = name.split('.')
        time_point = seq_info[3]
        time_point = time_point.split('-')[1]

        if time_point not in fasta_dict.keys():
            fasta_dict[time_point] = [fasta]
        else:
            fasta_dict[time_point] += [fasta]

    sequence_dict = {}

    #Put each time point into an array where the columns are sites and the rows
    #are sequences
    for time_point in fasta_dict.keys():
        curr_seq_list = fasta_dict[time_point]
        curr_array = []

        #loop through the seq.io objects
        for seq in curr_seq_list:
            curr_seq = seq.seq
            curr_array += [list(curr_seq)]
        
        #Convert the array to a numpy array
        curr_array = np.array(curr_array)
        sequence_dict[time_point] = curr_array

    #Make a dataframe to store the D' values
    d_prime_df = []


    #Calculate all of the D' values for each timepoint
    for time_point in sequence_dict.keys():
        curr_array = sequence_dict[time_point]
        seg_list = []
        

        #Transpose the array and loop through the sites
        curr_array = np.transpose(curr_array)
        for i in range(curr_array.shape[0]):
            curr_site = curr_array[i]
            curr_site = np.delete(curr_site, np.where(curr_site == '-'))

            #check if the site is segregating
            if len(set(curr_site)) > 1:
                seg_list += [i]
        
        #Calculate D' for each pair of segregating sites
        for i in range(len(seg_list)):
            for j in range(i):
                site_1 = seg_list[j]
                site_2 = seg_list[i]
                site_1_alleles = curr_array[site_1]
                site_2_alleles = curr_array[site_2]

                # print('**************************')
                # print(site_1_alleles)
                # print(site_2_alleles)

                #Remove gaps but preserve the haplotype information
                seq_gaps_1 = np.where(site_1_alleles == '-')
                site_1_alleles = np.delete(site_1_alleles, seq_gaps_1)
                site_2_alleles = np.delete(site_2_alleles, seq_gaps_1)

                seq_gaps_2 = np.where(site_2_alleles == '-')
                site_2_alleles = np.delete(site_2_alleles, seq_gaps_2)
                site_1_alleles = np.delete(site_1_alleles, seq_gaps_2)


                # print('------------------------')
                # print(site_1_alleles)
                # print(site_2_alleles)

                #Get the two most frequent alleles at site 1
                site_1_unique, site_1_counts = np.unique(site_1_alleles, return_counts=True)

                if len(site_1_unique) < 2:
                    continue

                #First save the most common then delete it and fine the second most common
                most_common = np.argmax(site_1_counts)
                allele_A = site_1_unique[most_common]
                site_1_unique = np.delete(site_1_unique, most_common)  
                site_1_counts = np.delete(site_1_counts, most_common)
                most_common = np.argmax(site_1_counts)
                allele_a = site_1_unique[most_common]

                #Get the two most frequent alleles at site 2
                site_2_unique, site_2_counts = np.unique(site_2_alleles, return_counts=True)

                if len(site_2_unique) < 2:
                    continue

                most_common = np.argmax(site_2_counts)
                allele_B = site_2_unique[most_common]
                site_2_unique = np.delete(site_2_unique, most_common)
                site_2_counts = np.delete(site_2_counts, most_common)
                most_common = np.argmax(site_2_counts)
                allele_b = site_2_unique[most_common]
                # print(allele_a, allele_A, allele_b, allele_B)

                #calculate the frequencies of the two most frequent alleles
                allele_A_freq = np.sum(site_1_alleles == allele_A) / len(site_1_alleles)
                allele_B_freq = np.sum(site_2_alleles == allele_B) / len(site_2_alleles)


                #Get the haplotype counts
                curr_haps = np.column_stack((site_1_alleles, site_2_alleles))
                curr_haps_unique, curr_haps_counts = np.unique(curr_haps, axis = 0, return_counts=True)

                #AB
                if [allele_A, allele_B] not in curr_haps_unique.tolist():
                    AB_count = 0
                else:
                    AB_index = np.where((curr_haps_unique[:,0] == allele_A) & (curr_haps_unique[:,1] == allele_B))[0][0]
                    AB_count = curr_haps_counts[AB_index]

                #Ab
                if [allele_A, allele_b] not in curr_haps_unique.tolist():
                    Ab_count = 0
                else:
                    Ab_index = np.where((curr_haps_unique[:,0] == allele_A) & (curr_haps_unique[:,1] == allele_b))[0][0]
                    Ab_count = curr_haps_counts[Ab_index]
                
                #aB
                if [allele_a, allele_B] not in curr_haps_unique.tolist():
                    aB_count = 0
                else:
                    aB_index = np.where((curr_haps_unique[:,0] == allele_a) & (curr_haps_unique[:,1] == allele_B))[0][0]
                    aB_count = curr_haps_counts[aB_index]
                
                #ab
                if [allele_a, allele_b] not in curr_haps_unique.tolist():
                    ab_count = 0
                else:
                    ab_index = np.where((curr_haps_unique[:,0] == allele_a) & (curr_haps_unique[:,1] == allele_b))[0][0]
                    ab_count = curr_haps_counts[ab_index]


                #put all of the D' values in a dataframe
                r_squared = r2.r2(AB_count, Ab_count, aB_count, ab_count)
                d_prime = r2.calc_D_prime(AB_count, Ab_count, aB_count, ab_count)
                d_val = r2.calc_D(AB_count, Ab_count, aB_count, ab_count)
                curr_time = timepoint_dict[time_point]
                
                d_prime_df.append([site_1, site_2, allele_A_freq, allele_B_freq,
                            AB_count, Ab_count, aB_count, ab_count,
                            r_squared, d_prime, d_val, curr_time])

    d_prime_df = pd.DataFrame(d_prime_df, columns = ['Locus_1', 'Locus_2', 'p_A', 'p_B',
                    'AB_obs', 'Ab_obs', 'aB_obs', 'ab_obs', 'r_squared', 'd_prime', 'd_val', 'timepoint'])        

    d_prime_df = d_prime_df[d_prime_df['p_A'].between(SEG_CUTOFF, 1-SEG_CUTOFF)]
    d_prime_df = d_prime_df[d_prime_df['p_B'].between(SEG_CUTOFF, 1-SEG_CUTOFF)]
    
    d_prime_df = d_prime_df[d_prime_df['timepoint'] > 0]

    if not os.path.exists(outDir + curr_dir + '/linkage/'):
        os.makedirs(outDir + curr_dir + '/linkage/')
    d_prime_df.to_pickle(outDir + curr_dir + '/linkage/r2_and_D')