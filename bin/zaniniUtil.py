import pandas as pd
import numpy as np
import sys
import os
import csv

def find_segregating_diagonal(coCounts_arr, all_seg = False):
    """ Takes a numpy array of co-SNP counts from the Zanini et al. data.
    Scans the diagonal of the array to find segregating sites. 
    ---------------------------------------------------------------------------
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

def getPatientList(dir):
    """ Takes a directory of files with patients and time points. 
    Generates a dictionary where the keys are patients and the values
    are lists of the corresponding data files.
    """
    participant_files = {}
    
    #loop through all of the vcf files and group them by participant
    for file in os.listdir(dir):
        #get the participant and the date
        parDate = file.split('.')[0]
        #get the participant
        curr_participant = parDate.split('_')[0]

        #if we have seen files for this participant already add them to entry
        if curr_participant in participant_files.keys():
            participant_files[curr_participant] += [parDate]
        else:
            participant_files[curr_participant] = [parDate]
    return participant_files

def label_vl_drats(stat_df, vlDir, cd4_label = False):
    """Labels a dataframe of D' ratios each with the average viral load between
    the two timepoints.
    ---------------------------------------------------------------------------
    Params
    ------------
    stat_df :   pd.DataFrame, containing the d' ratio for a pair of loci
                also includes d' at the first time point, and information 
                about the two timepoints of sampling
    vlDir :     str, a string containing the path to a directory with the viral
                load data for each individual in the form of tsv files
    cd4_label : bool, if True then the dataframe will be labeled with the CD4
                counts rather than the viral load measurmeents

    Returns
    -------
    labeled_rats: pd.DataFrame, input dataframe but with the viral loads 
                labeled
    """
    labeled_rats = []
    #label each participants ratios with the corresponding viral loads
    for curr_file in os.listdir(vlDir):
        dict_from_csv = {}
        if curr_file[0] == '.':
            continue

        #make a dictionary of timepoints and their viral loads
        with open(vlDir + curr_file, mode='r') as inp:
            reader = csv.reader(inp, delimiter= '\t')
            next(reader)
            dict_from_csv = {float(rows[0]):float(rows[1]) for rows in reader}
  
        #get the ratios for the participant
        participant = curr_file.split('.')[0]
        participant = participant.split('_')[1]
        par_d_rats = stat_df[stat_df['Participant'] == participant]

        #make a deep copy so we can set values on it
        curr_d_rats = par_d_rats.copy()
        curr_d_rats['Day_1'] = curr_d_rats['Time_1'] * 2
        curr_d_rats['Day_2'] = curr_d_rats['Time_2'] * 2
        

        #check if there are rows without any exact matrhing timepoints
        unmatched = curr_d_rats[~curr_d_rats['Day_1'].isin(dict_from_csv.keys()) | ~curr_d_rats['Day_2'].isin(dict_from_csv.keys())]
        if len(unmatched) > 0:
            print("Unmatched timepoints for participant " + participant)
            dict_from_csv = fix_unmatched_VLs_CD4(curr_d_rats, dict_from_csv)

        #label the ratios
        curr_d_rats = curr_d_rats[curr_d_rats['Day_1'].isin(dict_from_csv.keys())]
        curr_d_rats = curr_d_rats[curr_d_rats['Day_2'].isin(dict_from_csv.keys())]

        if cd4_label:
            curr_d_rats['CD4_1'] = curr_d_rats['Day_1'].map(lambda x: dict_from_csv[x])
            curr_d_rats['CD4_2'] = curr_d_rats['Day_2'].map(lambda x: dict_from_csv[x])
        else:
            curr_d_rats['VL_1'] = curr_d_rats['Day_1'].map(lambda x: dict_from_csv[x])
            curr_d_rats['VL_2'] = curr_d_rats['Day_2'].map(lambda x: dict_from_csv[x])

        #Check for viral load measurements between the two timepoints
        curr_d_rats = check_VL_between(curr_d_rats, dict_from_csv, cd4_label)
        labeled_rats.append(curr_d_rats)

    labeled_rats = pd.concat(labeled_rats, ignore_index= True)
    if cd4_label:
        labeled_rats['Ave_CD4'] = labeled_rats[['CD4_1', 'CD4_2']].mean(axis=1)
    else:
        labeled_rats['Ave_VL'] = labeled_rats[['VL_1', 'VL_2']].mean(axis=1)
    return labeled_rats

def combine_drats(d_rat_dir):
    """Combines D' ratio dataframes across individuals and fragments.
    ---------------------------------------------------------------------------
    Params
    ------------
    d_rat_dir: str, the path to the directory with the results of the snakemake
                estimation pipeline. 
    """
    stat_df = []
    for curr_data in os.listdir(d_rat_dir):
        if curr_data[0] == '.':
            continue
        run_info = curr_data.split('_')
        curr_par = run_info[0]
        curr_frag = run_info[1]

        d_ratio_file = d_rat_dir + curr_data + "/linkage/d_ratio"

        curr_stat_df = pd.read_pickle(d_ratio_file)
        curr_stat_df['Participant'] = curr_par
        curr_stat_df['Fragment'] = curr_frag
        stat_df.append(curr_stat_df)

    #put all the ratios together
    stat_df = pd.concat(stat_df, ignore_index= True)
    stat_df['d_i'] = stat_df['d_i'].to_numpy().astype(float)
    stat_df['d_ratio'] = stat_df['d_ratio'].to_numpy().astype(float)
    stat_df['Dist_X_Time'] = stat_df['Dist_X_Time'].to_numpy().astype(float)
    stat_df['d_i_1'] = stat_df['d_i'] * np.exp(np.negative(stat_df['d_ratio']))
    stat_df['Dist'] = stat_df['Locus_2'] - stat_df['Locus_1']
    return stat_df

########################### Helper Functions ##################################
def fix_unmatched_VLs_CD4(stat_df, label_dict):
    """For participants with no exact matching timepoints, pairs the closest 
    timepoints and labels the d' ratio.
    ---------------------------------------------------------------------------
    Params
    ------------
    stat_df :   pd.DataFrame, containing the d' ratio for a pair of loci
                also includes d' at the first time point, and information 
                about the two timepoints of sampling
    label_dict :   dict, a dictionary where the keys are the measurement 
                timepoints and the values are the viral loads at those times.
    Returns
    -------
    label_dict: dict, a dictionary where the keys are the measurement 
                timepoints and the values are the viral loads at those times.
                This dictionary has been augmented to include the VL of the 
                timepoint closest to the timepoint of interest for timepoints
                without an exact match.
    """
    #make a dataframe with the timepoints and viral loads
    keylist = label_dict.keys()

    for curr_time in stat_df['Day_1'].unique():
        #get the nearest timepoint
        closest_time = min(keylist, key=lambda x:abs(x-curr_time))
        #make sure the nearest timepoint is within 100 days
        if abs(closest_time - curr_time) < 100:
            label_dict[curr_time] = label_dict[closest_time]
    return label_dict


def check_VL_between(stat_df, label_dict, cd4_label):
    """Checks if there are any viral load samples between timepoints.
    ---------------------------------------------------------------------------
    Params
    ------------
    stat_df :   pd.DataFrame, containing the d' ratio for a pair of loci
                also includes d' at the first time point, and information 
                about the two timepoints of sampling
    label_dict :   dict, a dictionary where the keys are the measurement 
                timepoints and the values are the viral loads at those times.
    cd4_label : bool, if True then the dataframe will be labeled with the CD4
                counts rather than the viral load measurmeents
    Returns
    -------
    stat_df: pd.DataFrame, the same dataframe as the input but labeled with 
                any intermediate viral loads.
    """
    #make a dataframe with the timepoints and viral loads
    keylist = label_dict.keys()
    labelList = [label_dict[x] for x in keylist]
    if cd4_label:
        label_df = pd.DataFrame(zip(keylist, labelList), columns = ['Timepoint', 'CD4'])
        curr_label = 'CD4'
    else:
        label_df = pd.DataFrame(zip(keylist, labelList), columns = ['Timepoint', 'VL'])
        curr_label = 'VL'

    between_vls = []
    consec = []
    monoton = []
    for index, row in stat_df.iterrows():
        #look for timepoints in between the comparison timepoints
        row_btwn = label_df[label_df['Timepoint'].between(row['Day_1'], row['Day_2'], inclusive='neither')]

        #if there are timepoints between check if its monotonic and label the times between
        if len(row_btwn) > 0:
            #record that it's not consecutive and label the viral loads between
            curr_btwns = list(row_btwn[curr_label])
            between_vls.append(curr_btwns)
            consec.append(False)

            curr_btwns.extend([row[curr_label + '_1'], row[curr_label + '_2']])
            curr_vls = {row[curr_label + '_1'], row[curr_label + '_2']}
            if min(curr_btwns) in curr_vls \
                and max(curr_btwns) in curr_vls:
                    monoton.append(True)
            else: monoton.append(False)

        #otherwise just label everything as true automatically
        else:
            between_vls.append([])
            consec.append(True)
            monoton.append(True)
    stat_df['Between_VLs'] = between_vls
    stat_df['Consecutive'] = consec
    stat_df['Monotonic'] = monoton
    return stat_df

