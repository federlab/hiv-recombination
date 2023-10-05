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
################################ Neutral Data #################################
###############################################################################

neutral_data_dir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/revision/simulated_estimates_neutral/'

NUM_BOOTSTRAPS = 1000
NUM_GROUPS = 20
LOCI_PER_GROUP = 1000

all_stat_dfs_neut = pd.read_pickle(neutral_data_dir+ "all_stat_dfs_"+ \
                                        str(NUM_BOOTSTRAPS) + "_" + \
                                            str(NUM_GROUPS) + ".pkl")

LOCI_COUNTS = [] 
LDM_COUNTS = []

#loop through each rho value
for curr_rho in all_stat_dfs_neut['Sim_Rho'].unique():
    print(curr_rho)

    curr_rho_stat = all_stat_dfs_neut[all_stat_dfs_neut['Sim_Rho'] == curr_rho]

    #loop through each iter group
    for curr_iteration in range(0, NUM_GROUPS):
        #get the data for the current rho and iteration
        curr_stat_df = curr_rho_stat[
            curr_rho_stat['iter_group'] == curr_iteration]

        #downsample the group to a set number of segregating loci
        curr_stat_df = plne.downsample_loci(curr_stat_df, 
                                                LOCI_PER_GROUP)
        
        #count the number of segregating loci
        curr_count = slim_util.count_seg_sites(curr_stat_df)
        curr_count = curr_count['num_loci'].tolist()[0]
        LOCI_COUNTS.append([curr_count, curr_rho, curr_iteration])
        LDM_COUNTS.append([len(curr_stat_df), curr_rho, curr_iteration])


LOCI_COUNTS = pd.DataFrame(LOCI_COUNTS, columns = ['Segregating Sites', 
                                                'Sim_Rho', 'iter_group'])
LDM_COUNTS = pd.DataFrame(LDM_COUNTS, columns = ['LDMs',
                                                'Sim_Rho', 'iter_group'])

neut_mean_loc = round(LOCI_COUNTS['Segregating Sites'].mean())
neut_std_loc = round(LOCI_COUNTS['Segregating Sites'].std())
neut_mean_loc = f'{neut_mean_loc:,}' + r' pm ' + f'{neut_std_loc:,}'

neut_mean_ldm = round(LDM_COUNTS['LDMs'].mean())
neut_std_ldm = round(LDM_COUNTS['LDMs'].std())
neut_mean_ldm = f'{neut_mean_ldm:,}' + r' pm ' + f'{neut_std_ldm:,}'

my_table_df.append(['Neutral', neut_mean_loc, neut_mean_ldm])

###############################################################################
################################ Selection Data ###############################
###############################################################################

selection_data_dir = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/revision/simulated_estimates_selection/'

NUM_BOOTSTRAPS_sel = 1000
NUM_GROUPS_sel = 10

all_stat_dfs_sel = pd.read_pickle(selection_data_dir + "all_stat_dfs_"+ \
                                        str(NUM_BOOTSTRAPS_sel) + "_" + \
                                            str(NUM_GROUPS_sel) + ".pkl")

sel_loci_counts = slim_util.count_seg_sites(all_stat_dfs_sel)
sel_mean_loc = round(sel_loci_counts['num_loci'].mean())
sel_std_loc = round(sel_loci_counts['num_loci'].std())
sel_mean_loc = f'{sel_mean_loc:,}' + r' pm ' + f'{sel_std_loc:,}'

sel_mean_ldm = round(sel_loci_counts['num_ldms'].mean())
sel_std_ldm = round(sel_loci_counts['num_ldms'].std())
sel_mean_ldm = f'{sel_mean_ldm:,}' + r' pm ' + f'{sel_std_ldm:,}'


print(sel_loci_counts)
my_table_df.append(['Selection', sel_mean_loc, sel_mean_ldm])


###############################################################################
################################# In Vivo Data ################################
###############################################################################
inv_data_dir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/"
vlDir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/"

#Make the dataframe containing D' ratios
og_stat_df = zu.combine_drats(inv_data_dir )
og_stat_df = zu.label_vl_drats(og_stat_df, vlDir)

og_stat_df = og_stat_df[og_stat_df['Fragment'] != 'F5']
og_stat_df = og_stat_df[og_stat_df['Time_Diff'].gt(50)]
QUANTILE_LIST = [0, 0.33, 0.66, 1]

lower_label = "Lower Tertile (Fig. 5 Data)"
middle_label = "Middle Tertile (Fig. 5 Data)" 
upper_label = "Upper Tertile (Fig. 5 Data)"

#Only get values for which the viral load doesn't leave the boundaries of the 
#initial and final measurments
stat_df = og_stat_df[og_stat_df['Monotonic'] == True]

############################### Tertiles All ##################################
stat_df['VL_Tertile'], bins = pd.qcut(og_stat_df['Ave_VL'], 
                            q = QUANTILE_LIST, retbins= True, 
                            labels = [lower_label, middle_label, upper_label])
grouped_stats = stat_df.groupby(by = ['VL_Tertile'])

#Loop through each of the groups and count the segregating loci
for name, group in grouped_stats:
    curr_count = zu.count_seg_sites(group, 'VL_Tertile')
    
    #Get the mean and standard deviation of the number of segregating loci
    curr_loc = round(curr_count['num_loci'].tolist()[0])
    curr_loc = f'{curr_loc:,}'

    #Get the mean and standard deviation of the number of LDMs
    curr_ldm = round(curr_count['num_ldms'].tolist()[0])
    curr_ldm = f'{curr_ldm:,}'
    
    my_table_df.append([name, curr_loc, curr_ldm])

################################### Half P1 ###################################
p1_stat_df = og_stat_df[og_stat_df['Participant'] == 'p1']


my_labels = ["Upper Half (Fig. 6 Data P1)","Lower Half (Fig. 6 Data P1)"]
p1_stat_df['VL_Quantile'], bins = pd.qcut(p1_stat_df['Ave_VL'], 
                                q = [0, 0.5, 1], retbins= True,
                                labels = my_labels)

grouped_stats = p1_stat_df.groupby(by = ['VL_Quantile'])

#Loop through each of the groups and count the segregating loci
for name, group in grouped_stats:
    curr_count = zu.count_seg_sites(group, 'VL_Quantile')
    
    #Get the mean and standard deviation of the number of segregating loci
    curr_loc = round(curr_count['num_loci'].tolist()[0])
    curr_loc = f'{curr_loc:,}'

    #Get the mean and standard deviation of the number of LDMs
    curr_ldm = round(curr_count['num_ldms'].tolist()[0])
    curr_ldm = f'{curr_ldm:,}'
    
    my_table_df.append([name, curr_loc, curr_ldm])


################################### Half P3 ###################################
p3_stat_df = og_stat_df[og_stat_df['Participant'] == 'p3']


my_labels = ["Upper Half (Fig. 6 Data P3)","Lower Half (Fig. 6 Data P3)"]
p3_stat_df['VL_Quantile'], bins = pd.qcut(p3_stat_df['Ave_VL'], 
                                q = [0, 0.5, 1], retbins= True,
                                labels = my_labels)

grouped_stats = p3_stat_df.groupby(by = ['VL_Quantile'])

#Loop through each of the groups and count the segregating loci
for name, group in grouped_stats:
    curr_count = zu.count_seg_sites(group, 'VL_Quantile')
    
    #Get the mean and standard deviation of the number of segregating loci
    curr_loc = round(curr_count['num_loci'].tolist()[0])
    curr_loc = f'{curr_loc:,}'

    #Get the mean and standard deviation of the number of LDMs
    curr_ldm = round(curr_count['num_ldms'].tolist()[0])
    curr_ldm = f'{curr_ldm:,}'
    
    my_table_df.append([name, curr_loc, curr_ldm])

################################### Half P4 ###################################
p4_stat_df = og_stat_df[og_stat_df['Participant'] == 'p4']


my_labels = ["Upper Half (Fig. 6 Data P4)","Lower Half (Fig. 6 Data P4)"]
p4_stat_df['VL_Quantile'], bins = pd.qcut(p4_stat_df['Ave_VL'], 
                                q = [0, 0.5, 1], retbins= True,
                                labels = my_labels)

grouped_stats = p4_stat_df.groupby(by = ['VL_Quantile'])

#Loop through each of the groups and count the segregating loci
for name, group in grouped_stats:
    curr_count = zu.count_seg_sites(group, 'VL_Quantile')
    
    #Get the mean and standard deviation of the number of segregating loci
    curr_loc = round(curr_count['num_loci'].tolist()[0])
    curr_loc = f'{curr_loc:,}'

    #Get the mean and standard deviation of the number of LDMs
    curr_ldm = round(curr_count['num_ldms'].tolist()[0])
    curr_ldm = f'{curr_ldm:,}'
    
    my_table_df.append([name, curr_loc, curr_ldm])

################################### Half P7 ###################################
p7_stat_df = og_stat_df[og_stat_df['Participant'] == 'p7']


my_labels = ["Upper Half (Fig. 6 Data P7)","Lower Half (Fig 6. Data P7)"]
p7_stat_df['VL_Quantile'], bins = pd.qcut(p7_stat_df['Ave_VL'], 
                                q = [0, 0.5, 1], retbins= True,
                                labels = my_labels)

grouped_stats = p7_stat_df.groupby(by = ['VL_Quantile'])

#Loop through each of the groups and count the segregating loci
for name, group in grouped_stats:
    curr_count = zu.count_seg_sites(group, 'VL_Quantile')
    
    #Get the mean and standard deviation of the number of segregating loci
    curr_loc = round(curr_count['num_loci'].tolist()[0])
    curr_loc = f'{curr_loc:,}'

    #Get the mean and standard deviation of the number of LDMs
    curr_ldm = round(curr_count['num_ldms'].tolist()[0])
    curr_ldm = f'{curr_ldm:,}'
    
    my_table_df.append([name, curr_loc, curr_ldm])

my_table_df = pd.DataFrame(my_table_df, columns = ['Estimation Group',
                                        r'# of Segregating Sites',
                                        r'# of LDMs'])

print(my_table_df.to_latex(index=False))

