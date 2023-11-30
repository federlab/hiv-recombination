import os
import sys
import numpy as np
import pandas as pd

#This file contains helper functions for working with the simulated data

############################### Labeling Functions ############################
def label_rho(my_dataframe):
    """ This function takes a dataframe and labels the rho values with floats
    instead of strings so that they can be used as a continuous variable
    To get the comparison, we need to convert the rho values to floats."""
    rhoDict = {"0.001" : 0.001,
               "0.002" : 0.002,
               "0.004" : 0.004,
               "0.005" : 0.005,
                "1e-04" : 0.0001,
                "0.00011": 0.00011,
                "0.00012": 0.00012,
                "0.00015": 0.00015,
                "2e-04" : 0.0002,
                "4e-04" : 0.0004,
                "5e-04" : 0.0005,
                "1e-05" : 0.00001,
                "1.1e-05" : 0.000011,
                "1.2e-05" : 0.000012,
                "1.5e-05" : 0.000015,
                "2e-05" : 0.00002,
                "4e-05" : 0.00004,
                "5e-05" : 0.00005,
                "1e-06" : 0.000001,
                "1.1e-06" : 0.0000011,
                "1.2e-06" : 0.0000012,
                "1.5e-06" : 0.0000015,
                "2e-06" : 0.000002,
                "4e-06" : 0.000004,
                "5e-06" : 0.000005}
    rho_dict_fix_strings = { "0.001" : r"$10^{-3}$",
                            "0.002" : r"$2\times10^{-3}$",
                            "0.004" : r"$4\times10^{-3}$",
                            "0.005" : r"$5\times10^{-3}$",
                            "1e-04" : r"$10^{-4}$",
                            "0.00011": r"$1.1\times10^{-4}$",
                            "0.00012": r"$1.2\times10^{-4}$",
                            "0.00015": r"$1.5\times10^{-4}$",
                            "2e-04" : r"$2\times10^{-4}$",
                            "4e-04" : r"$4\times10^{-4}$",
                            "5e-04" : r"$5\times10^{-4}$",
                            "1e-05" : r"$10^{-5}$",
                            "1.1e-05" : r"$1.1\times10^{-5}$",
                            "1.2e-05" : r"$1.2\times10^{-5}$",
                            "1.5e-05" : r"$1.5\times10^{-5}$",
                            "2e-05" : r"$2\times10^{-5}$",
                            "4e-05" : r"$4\times10^{-5}$",
                            "5e-05" : r"$5\times10^{-5}$",
                            "1e-06" : r"$10^{-6}$",
                            "1.1e-06" : r"$1.1\times10^{-6}$",
                            "1.2e-06" : r"$1.2\times10^{-6}$",
                            "1.5e-06" : r"$1.5\times10^{-6}$",
                            "2e-06" : r"$2\times10^{-6}$",
                            "4e-06" : r"$4\times10^{-6}$",
                            "5e-06" : r"$5\times10^{-6}$"}


    #redo the labeling on the rho values from what was used in the simulation
    intRhoList = []
    newStringRho = []
    for entry in my_dataframe['Sim_Rho']:
        intRhoList.append(rhoDict[entry])
        newStringRho.append(rho_dict_fix_strings[entry])
    my_dataframe['Sim_float_rho'] = intRhoList
    my_dataframe['Sim_Rho'] = newStringRho

    return my_dataframe

def label_m_rate(my_dataframe):
    """ This function takes in a dataframe with a migration rate column titled
    'm_rate'. Then, it converts the string to scientific notation for better
    labeling and returns the newly labeled dataframe
    """
    m_dict_fix_strings = {
            "0.01" : r"$m=10^{-2}$",
            "0.001" : r"$m=10^{-3}$",
            "1e-04" : r"$m=10^{-4}$",
            "1e-05" : r"$m=10^{-5}$",
            "Well Mixed" : "Well Mixed"}
    newStringM = []
    for entry in my_dataframe['m_rate']:
        newStringM.append(m_dict_fix_strings[entry])
    my_dataframe['m_rate'] = newStringM

    return my_dataframe

def label_Ne(my_dataframe):
    """ This function takes in a dataframe with a Ne column titled
    'Ne'. Then, it converts the string to scientific notation for better
    labeling and returns the newly labeled dataframe
    """
    ne_dict_fix_strings = {
            "1000" : r"$N_e=10^3$",
            "5000" : r"$N_e=5\times10^3$",
            "10000" : r"$N_e=10^4$",
            "50000" : r"$N_e=5\times10^4$",
            "1e+05" : r"$N_e=10^5$",
            "5e+05" : r"$N_e=5\times10^5$"}
    newStringNe = []
    for entry in my_dataframe['Ne']:
        newStringNe.append(ne_dict_fix_strings[entry])
    my_dataframe['Ne'] = newStringNe

    return my_dataframe


################## Functions for calculating statistics #######################
def calc_pop_metrics(slim_out_file, genome_length):
    """ This function takes in a slim output file and returns the nucleotide
    diversity (mean pairwise hamming distance) and mean sequence divergence for
    plotting.
    ---------------------------------------------------------------------------
    params:
    ---------
    slim_out_file:  str, path to the slim output file we'll be calculating the 
                        metrics for
    genome_length:  int, length of the genome that was simulated (in bp)
    
    returns:
    ---------
    dist_df:    pd.DataFrame, dataframe containing the mean pairwise distance,
                timepoint, and mean sequence divergence
    """
    #Make a list where we'll store the results
    dist_df = []

    #Read in the slim output
    with open(slim_out_file) as f:
        lines = f.read()
        #Separate each output section and get rid of the script text
        lines = lines.split('#OUT')[1:]

    #Now separately parse each section
    for curr_section in lines:
        #Get the timepoint
        first_line = curr_section.split('\n', 1)[0]
        first_line = first_line.split(' ')

        #Isolate the timepoint and convert to an int
        time = int(first_line[1])

        #Check if there is this section includes genotype info
        contains_genomes = len(curr_section.split('Genomes:')) > 1
        
        #Separate out and process the genotype info if it exists
        if contains_genomes:
            #Separate out the genotype info
            split_results = curr_section.split('Genomes:')
            genome_data = split_results[1]

            #Make a list where each entry is the list of mutations in a genome
            genome_data = genome_data.split('\n')[1:]
            genome_data = [x.split(' ')[2:] for x in genome_data]

            pairwise_dists = []
            seq_div = []

            #Make all pairwise comparisons
            for j in range(len(genome_data)):
                seq_1 = genome_data[j]
                seq_div.append(len(seq_1)/genome_length)
                for k in range(j+1, len(genome_data)):
                    seq_2 = genome_data[k]

                    #Get the pairwise distance
                    total_muts = len(seq_1) + len(seq_2)
                    shared_muts = set(seq_1 + seq_2)
                    pairwise_dists.append(
                        (total_muts - 
                            (total_muts - len(shared_muts)))/genome_length)
            
            #Get the mean pairwise distance
            mean_pairwise_dist = np.mean(pairwise_dists)
            mean_seq_div = np.mean(seq_div)
            dist_df.append([time, mean_pairwise_dist, mean_seq_div])

    dist_df = pd.DataFrame(dist_df, 
            columns = ['timepoint', 'mean_pairwise_dist', 'mean_sequence_div']) 

    return dist_df

def count_seg_sites(stat_df):
    """ This function takes the D' ratio dataframe used to produce grouped
    estimates of rho. Then, it calculates the number of segregating sites that
    went into each estimate.
    ---------------------------------------------------------------------------
    params:
    ---------
    stat_df:   pd.DataFrame, dataframe containing the D' ratio values output by
                the estimation_util module

    returns:
    ---------
    all_loci_counts:    pd.DataFrame, dataframe containing the number of loci
                        for each estimate group
    """
    all_loci_counts = []

    #Loop through each rho value
    for curr_rho in stat_df['Sim_Rho'].unique():
        curr_rho_stat = stat_df[stat_df['Sim_Rho'] == curr_rho]

        #Loop through each group we were estimating from
        for curr_group in curr_rho_stat['iter_group'].unique():

            #get the data for the current rho and iteration
            curr_stat_df = curr_rho_stat[
                    curr_rho_stat['iter_group'] == curr_group]

            #make a place to store the number of loci in the estimate
            group_num_loci = 0

            #Next we need to count segregating sites but we have to group by
            # locus and simulation
            for curr_rep in curr_stat_df['rep'].unique():
                rep_stat_df = curr_stat_df[curr_stat_df['rep'] == curr_rep]
                curr_loci = set(rep_stat_df['Locus_1'].unique().tolist() + \
                                rep_stat_df['Locus_2'].unique().tolist())
                group_num_loci += len(curr_loci)
            num_ldms = len(curr_stat_df)

            #Now we have the number of loci in the group, so we can add it to 
            # the dataframe
            all_loci_counts.append([curr_rho, curr_group, group_num_loci,
                                    num_ldms])

    #Make a dataframe of the counts
    all_loci_counts = pd.DataFrame(all_loci_counts, 
                columns = ['Sim_Rho', 'iter_group', 'num_loci', 'num_ldms'])
    
    return all_loci_counts

def make_dprime_df(dataDir, allele_freq_threshes, distance_thresh,
                   four_hap_thresh = True, migration = True, rho_list = None,
                   NE_VAL = '10000', rep_limit = 250):
    """This function makes a dataframe of the D' ratios labeled by rho value 
    and migration rate. I am going to use it to make linkage plots 
    particularly for the population substructure data. 
    ---------------------------------------------------------------------------
    params:
    ---------
    dataDir:    str, path to the directory containing the slim output files
    allele_freq_thresh: tuple, the minimum and maximum allele frequencies that
                we'll use to calculate D' ratios
    distance_thresh: int, the maximum distance between loci that we'll put in
                the dataframe
    four_hap_thresh: bool, whether to require four haplotypes to include ratios
                in the dataframe
    migration:  bool, whether to include the migration rate in the dataframe
    rho_list:  list, list of rho values to include in the dataframe
    NE_VAL:     str, the Ne value for the data to include in the dataframe
    rep_limit:  int, the maximum number of replicates to include in the
                dataframe

    returns:
    ---------
    dprime_df:  pd.DataFrame, dataframe containing the D' ratios for all of 
                the simulations
    """
    #Make a list where we'll store the results
    dprime_df = []

    #Loop through each dataset in the directory
    for curr_dir in os.listdir(dataDir):
        #Skip the non-directories
        if not os.path.isdir(dataDir+curr_dir) or curr_dir[0] == '.':
            continue
        
        # print(curr_dir, file = sys.stderr)
        #Get the info for the run
        run_info = curr_dir.split('_')
        rho_val = run_info[1]
        rho_val = rho_val[3:]
        ne_val = run_info[2]
        ne_val = ne_val[2:]
        rep = run_info[4]
        rep = rep[3:]

        #Check if we want to include this data based on it's rho value
        if rho_list:
            if rho_val not in rho_list:
                continue
        if ne_val != NE_VAL:
            continue

        if int(rep) > rep_limit-1:
            continue


        #Now, get the D' ratios for the run
        curr_file = dataDir + curr_dir + "/linkage/r2_and_D"
        D_vals_df = pd.read_pickle(curr_file)
        D_vals_df['Rho'] = rho_val
        D_vals_df['Rep'] = rep        

        #Get the migration rate if it exists
        if migration:
            m_rate = run_info[5]
            m_rate = m_rate[2:]
            D_vals_df['m_rate'] = m_rate

        #Filter the data by allele frequency
        freq_min = allele_freq_threshes[0]
        freq_max = allele_freq_threshes[1]
        D_vals_df = D_vals_df[D_vals_df['p_A'].between(freq_min, freq_max)]
        D_vals_df = D_vals_df[D_vals_df['p_B'].between(freq_min, freq_max)]

        #Filter the data by distance
        D_vals_df['Dist'] = D_vals_df['Locus_2'] - D_vals_df['Locus_1']
        D_vals_df = D_vals_df[D_vals_df['Dist'] <= distance_thresh]

        #Filter the data by number of haplotypes
        if four_hap_thresh:
            D_vals_df = D_vals_df[D_vals_df['AB_obs'] > 0]
            D_vals_df = D_vals_df[D_vals_df['Ab_obs'] > 0]
            D_vals_df = D_vals_df[D_vals_df['aB_obs'] > 0]
            D_vals_df = D_vals_df[D_vals_df['ab_obs'] > 0]

        #Add the data to the list we'll return
        dprime_df.append(D_vals_df)
    
    dprime_df = pd.concat(dprime_df, ignore_index=True)

    return dprime_df

