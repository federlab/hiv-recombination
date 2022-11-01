import numpy as np
import pandas as pd


def downsample_genotype_df(genotype_df, num_genomes, dist_samples, min_read_len, max_read_len, half_read_len):
    """ Takes in a genotype dataframe of simulated data and downsamples
    the read coverage so that it matches the empirical distribution of
    read coverage based on the Zanini et al. paired end sampling
    ---------------------------------------------------------------------------
    params:
    ---------
    genotype_df:    pd.DataFrame, genotype dataframe of simulated data
    num_genomes:    int, number of genomes to downsample to
    dist_samples:   int, number of samples to use to make the empirical
                        distribution
    min_read_len:   int, minimum read length possible
    max_read_len:   int, maximum read length possible
    half_read_len:  int, half the read length

    returns:
    ---------
    downsampled_df: pd.DataFrame, downsampled genotype dataframe
    """
    #Make the empirical distribution of coverage
    coverage_dist = make_downsampled_read_dist(dist_samples, min_read_len, \
                        max_read_len, half_read_len)

    #Group the genotype dataframe by position
    genotype_df = genotype_df.groupby(['Locus_1', 'Locus_2', 'timepoint'])

    downsampled_df = []
    #Downsample each pair to the empirically expected coverage
    for name, group in genotype_df:
        distance = abs(name[0] - name[1])

        #Get the expected coverage
        if distance < 0 or distance > max_read_len:
            curr_coverage = 0
        else:
            curr_coverage = coverage_dist.iloc[[distance]]
            curr_coverage = curr_coverage['coverage'].tolist()[0]

        #Downsample the group to the expected coverage
        group_total = np.sum(group['Count'])
        group_probs = np.array(group['Count']) / group_total
        new_counts = np.random.multinomial(curr_coverage * num_genomes, group_probs)

        #Update the genotype dataframe
        new_count_df = group.copy()
        new_count_df['Count'] = new_counts
        downsampled_df.append(new_count_df)

    downsampled_df = pd.concat(downsampled_df, ignore_index=True)
    return downsampled_df        

###############################################################################
########################### Helper Functions ##################################

def make_downsampled_read_dist(dist_samples, min_read_len, max_read_len, half_read_len):
    """ Makes an empirical distribution of coverage based on the distance 
    between two loci and the paired-end sampling of reads.
    ---------------------------------------------------------------------------
    params:
    ---------
    dist_samples:    int, number of samples to use to make the empirical 
                        distribution
    min_read_len:   int, minimum read length possible
    max_read_len:   int, maximum read length possible
    half_read_len:  int, half the read length
    
    returns:
    ---------
    coverage_df:    pd.DataFrame, empirical distribution of coverage with 
                    the distance column indicating the distance between
                    loci and the coverage column indicating proportion of
                    coverage at that distance.
    """
    #First get the start positions of the reads by sampling uniformly across
    # positions
    read_starts = np.floor(np.random.uniform(0, max_read_len, dist_samples))

    #Draw insert sizes for each of the reads
    #we know the distribution of read lengths and that each half of a read 
    #is a set number of bp, but we don't know the distribution of insert sizes
    insert_sizes = np.floor(
                    np.random.uniform(min_read_len, max_read_len, dist_samples)\
                                - 2 * half_read_len)
    
    #Insert sizes can't be negative
    insert_sizes = np.where(insert_sizes < 0, 0, insert_sizes)

    #Make a dataframe of reads giving the positions of the two halves
    sim_reads = pd.DataFrame(list(zip(read_starts, 
                    read_starts + half_read_len,
                    read_starts + half_read_len + insert_sizes,
                    read_starts + 2 * half_read_len + insert_sizes)),
                    columns = ['L_start', 'L_end', 'R_start', 'R_end'])

    #Make sure reads all contain the central position (which is arbitrarily
    #the same as the max read length)
    sim_reads = sim_reads[\
        ((sim_reads['L_start'] <= max_read_len) & (sim_reads['L_end'] >= max_read_len)) \
        | ((sim_reads['R_start'] <= max_read_len) & (sim_reads['R_end'] >= max_read_len))]

    read_coverage_arr = np.zeros(2 * max_read_len + 1)
    #For all possible positions, check if it is in either portion of the read
    for i in range(2 * max_read_len + 1):
        #get all of the reads containing i
        contains_i = sim_reads[((sim_reads['L_start'] <= i) & (sim_reads['L_end'] >= i)) \
                        | ((sim_reads['R_start'] <= i) & (sim_reads['R_end'] >= i))]

        #divide by the total number of reads
        read_coverage_arr[i] = len(contains_i) / len(sim_reads)


    #make an array of all of the positions relative to the focal position
    position_arr = np.arange(2 * max_read_len + 1)
    position_arr = list(map(lambda x: abs(max_read_len - x), position_arr))
    
    coverage_df = pd.DataFrame(list(zip(position_arr, read_coverage_arr)),
                                columns = ['distance', 'coverage'])

    coverage_df = coverage_df.groupby('distance')
    coverage_df = coverage_df.aggregate(np.mean)
    
    return coverage_df

