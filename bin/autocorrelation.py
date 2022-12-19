import pandas as pd
import numpy as np
import r2Analysis as r2
import sys

#Estimate the recombination rate by looking at the decay of 
#linkage over time. We need to look at the autocorrelation of linkage over time
def calculate_d_ratios(linkage_file, threshold = 0.2, four_haps = True, stat = 'd_prime'):
    """ Takes in a dataframe of D' statistics. For each pair of timepoints and
    each pair of loci, calculates the ratio of D' values at the two timepoints.
    ---------------------------------------------------------------------------
    Params
    -------
    linkage_file :    str, path to the dataframe with D' statistics
    threshold :       float, D' ratios are only included in the output if the D'
                        value at the first timepoint is > THRESHOLD
    num_haps :        bool, the number of haplotypes required for each timepoint
                      to be included. If True, requires 4 haplotypes, if False
                      only requires 3 haplotypes
    stat :            str, the statistic to use for the calculation. Can be either
                     'd_prime' or 'd_val'
    Returns
    -------
    stat_df :   pd.DataFrame, containing the d' ratio for a pair of loci
                also includes d' at the first time point, and information 
                about the two timepoints of sampling
    """
    #first just try and print the array
    rd_arr = pd.read_pickle(linkage_file)
    rd_arr.dropna(inplace = True)
    rd_arr['support'] = rd_arr['AB_obs'] + rd_arr['Ab_obs'] + \
         rd_arr['aB_obs'] + rd_arr['ab_obs']
    
    #Only include loci with at least three haplotypes
    three_filter = []
    for i in range(len(rd_arr)):
        my_row = rd_arr.iloc[i]
        row_sum = 0
        if my_row['AB_obs'] > 0:
            row_sum += 1
        if my_row['aB_obs'] > 0:
            row_sum += 1
        if my_row['Ab_obs'] > 0:
            row_sum += 1
        if my_row['ab_obs'] > 0:
            row_sum += 1
        three_filter.append(row_sum)
    rd_arr['three_filter'] = three_filter

    print(((rd_arr['AB_obs'] > 0) + (rd_arr['aB_obs'] > 0) \
            + (rd_arr['Ab_obs'] > 0) + (rd_arr['ab_obs'] > 0)), file = sys.stderr)
    #only include loci with 4 haplotypes
    if four_haps:
        rd_arr = rd_arr.loc[(rd_arr['AB_obs'] > 0) & (rd_arr['aB_obs'] > 0) \
            & (rd_arr['Ab_obs'] > 0)  & (rd_arr['ab_obs'] > 0) ]
    else:
        rd_arr = rd_arr.loc[rd_arr['three_filter'] > 2 ]

    #group by the pair of loci
    grouped_loci = rd_arr.groupby(['Locus_1', 'Locus_2'])

    #make a dataframe to save the results in
    stat_df = []

    #loop over the pairs of loci to calculate -log(D_(t_i+1)/D_(t_i))
    for name, group in grouped_loci:
        #loop over the timepoints
        group_times = group['timepoint'].unique()
        group_times.sort()
        for i in range(len(group_times) - 1):

            curr_time = float(group_times[i])
            next_time = float(group_times[i+1])

            d_i = group[group['timepoint'] == curr_time]
            d_i = d_i[stat].tolist()[0]

            d_i_1 = group[group['timepoint'] == next_time]
            d_i_1 = (d_i_1[stat]).tolist()[0]
            
            if d_i < threshold:
                # print('D_i is zero')
                continue
            if d_i == 0:
                continue

            curr_val = -np.log(d_i_1/d_i)
            stat_df.append([curr_val, name[0], name[1], next_time - curr_time, d_i, curr_time, next_time])

    #Add necessary info to statistic dataframe and only get entries with relevant Dist_X_Time
    stat_df = pd.DataFrame(stat_df, columns = ['d_ratio', 'Locus_1', 'Locus_2', 'Time_Diff', 'd_i', 'Time_1', 'Time_2'])
    stat_df['Dist_X_Time'] = (stat_df['Locus_2'] - stat_df['Locus_1']) * stat_df['Time_Diff']
    stat_df = stat_df[stat_df['d_ratio'].between(-10,10)]
    stat_df.replace([np.inf, -np.inf], np.nan, inplace=True)
    stat_df.dropna(inplace = True)

    return stat_df
