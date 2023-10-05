import sys
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import numpy as np
import pandas as pd
import slimUtil as slim_util
import plot_neher as plne
import zaniniUtil as zu

#This table will record the numbers of segregating sites used for all of the estimates displayed in the paper

my_table_df = []

###############################################################################
################################# In Vivo Data ################################
###############################################################################
inv_data_dir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/"
vlDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/"

#Make the dataframe containing D' ratios
stat_df = zu.combine_drats(inv_data_dir )
stat_df = zu.label_vl_drats(stat_df, vlDir)

stat_df = stat_df[stat_df['Fragment'] != 'F5']
stat_df = stat_df[stat_df['Time_Diff'].gt(50)]
stat_df = stat_df[stat_df['Participant'] != 'p10']

QUANTILE_LIST = [0, 0.33, 0.66, 1]


#Only get values for which the viral load doesn't leave the boundaries of the 
#initial and final measurments
stat_df = stat_df[stat_df['Monotonic'] == True]

stat_df['Participant'] = stat_df['Participant'].str.replace('p','Participant ')

grouped_stats = stat_df.groupby(by = ['Participant'])

#Loop through each of the groups and count the segregating loci
for name, group in grouped_stats:
    curr_count = zu.count_seg_sites(group, 'Participant')
    
    #Get the mean and standard deviation of the number of segregating loci
    curr_loc = round(curr_count['num_loci'].tolist()[0])
    curr_loc = f'{curr_loc:,}'

    #Get the mean and standard deviation of the number of LDMs
    curr_ldm = round(curr_count['num_ldms'].tolist()[0])
    curr_ldm = f'{curr_ldm:,}'
    
    my_table_df.append([name, curr_loc, curr_ldm, int(name.split(' ')[-1])])

my_table_df = pd.DataFrame(my_table_df, columns = [' ',
                                        r'# of Segregating Sites',
                                        r'# of LDMs', 'sortkey'])
my_table_df = my_table_df.sort_values(by = ['sortkey'])
my_table_df = my_table_df.drop(columns = ['sortkey'])


print(my_table_df.to_latex(index=False))
