import os
import numpy as np
import pandas as pd
import plot_neher as plne

# This file can handle dividing the reps into groups and estimating recombination rates

def make_comparison_dataframes(dataDir, outDir, NUM_REPS, NUM_GROUPS,
                                NUM_BOOTSTRAPS, DIST_TIME_MAX, THRESHOLD = 0.2,
                                DSTATISTIC = False, DOWNSAMPLE = True,
                                LOCI_PER_GROUP = None, PREMADE_DF = None):
    """This function randomly pairs simulations and makes a dataframe
    that shows how often we can correctly order pairs and how often
    their confidence intervals overlap.
    ---------------------------------------------------------------------------
    Parameters:
    -----------
    dataDir:    str, the directory holding the data we will use for estimation.
    outDir:     str, the directory for saving the output dataframes.
    NUM_REPS:   int, the number of reps to divide into groups.
    NUM_GROUPS: int, the number of groups we will use for the analysis.
    NUM_BOOTSTRAPS: int, the number of bootstraps to use for the analysis.
    DIST_TIME_MAX:  int, the maximum distance * time to fit the curves up to.
    THRESHOLD:  float, the minimum D' value at timepoint 1 to include in the
                    analysis.
    DSTATISTIC: bool, whether or not to use the D statistic rather than D'
    DOWNSAMPLE: bool, whether or not to downsample the datasets.
    LOCI_PER_GROUP:   int, the number of segregating loci to downsample each
                        group to. Only used if DOWNSAMPLE is True.
    PREMADE_DF: str, filepath for existing d_stat dataframe to estimate from.

    Returns:
    ------------
    None, but saves the estimate and confidence interval dataframes to
         the outDir.
    """
    if not PREMADE_DF:
        #First I am going to read the data and randomly pair simulations
        all_stat_dfs = []
        
        #loop through each of the dataframes for the separate simulations
        for curr_data in os.listdir(dataDir):
            #only get the data directories, not hidden files
            if curr_data[0] == '.':
                continue

            #get the information for the current run
            run_info = curr_data.split('_')
            sim_rho = run_info[1]
            sim_rho = sim_rho[3:]
            rep = run_info[-1]

            if int(rep[3:]) >= NUM_REPS:
                continue


            #get the dataframe for the current run
            if DSTATISTIC:
                d_ratio_file = dataDir + curr_data + "/linkage_D/d_ratio_three_haps"
            else:
                d_ratio_file = dataDir + curr_data + "/linkage/d_ratio"

            stat_df = pd.read_pickle(d_ratio_file)

            stat_df['rep'] = int(rep[3:])
            stat_df['Sim_Rho'] = sim_rho
            all_stat_dfs.append(stat_df)
        all_stat_dfs = pd.concat(all_stat_dfs)

        #Filter the groups so the first d' value is greater than 0.2
        all_stat_dfs = all_stat_dfs[all_stat_dfs['d_i'] > THRESHOLD]
    else:
        all_stat_dfs = pd.read_pickle(PREMADE_DF)

    #Randomly divide the reps into 10 groups
    rep_groups = np.array(range(0, NUM_REPS))
    np.random.shuffle(rep_groups)
    rep_groups = np.array_split(rep_groups, NUM_GROUPS)


    #Make a dictionary to label each group
    group_dict = {}
    for i in range(len(rep_groups)):
        for j in rep_groups[i]:
            group_dict[j] = i
            
    group_labels = [group_dict[x] for x in all_stat_dfs['rep']]
    all_stat_dfs['iter_group'] = group_labels

    print("Loaded the data Successfully")

    ########################## Estimating recombination rates ##################
    #loop through each of the distance cutoffs
    all_stat_dfs = all_stat_dfs[all_stat_dfs['Dist_X_Time'] <= DIST_TIME_MAX]
    all_stat_dfs.to_pickle(outDir + "all_stat_dfs_" + \
                               str(NUM_BOOTSTRAPS) + "_" + str(NUM_GROUPS) \
                                + ".pkl")
    print(all_stat_dfs.head())
    all_estimate_df = []

    #loop through each rho value
    for curr_rho in all_stat_dfs['Sim_Rho'].unique():
        print(curr_rho)
        curr_rho_stat = all_stat_dfs[all_stat_dfs['Sim_Rho'] == curr_rho]

        #loop through each iter group
        for curr_iteration in range(0, NUM_GROUPS):
            #get the data for the current rho and iteration
            curr_stat_df = curr_rho_stat[
                curr_rho_stat['iter_group'] == curr_iteration]

            #downsample the group to a set number of segregating loci
            if DOWNSAMPLE:
                curr_stat_df = plne.downsample_loci(curr_stat_df, 
                                                    LOCI_PER_GROUP)

            lower_fit, upper_fit, estimate_df = \
                plne.bootstrap_rho(curr_stat_df, NUM_BOOTSTRAPS)

            estimate_df['Group'] = curr_iteration
            estimate_df['Sim_Rho'] = curr_rho
            all_estimate_df.append(estimate_df)



    all_estimate_df = pd.concat(all_estimate_df, ignore_index=True)
    all_estimate_df.to_pickle(outDir + "all_estimate_df_" + \
                               str(NUM_BOOTSTRAPS) + "_" + str(NUM_GROUPS)\
                                + ".pkl")

    all_conf_ints = []

    #Drop NA values before we get the confidence interval so they don't interfere
    all_estimate_df = all_estimate_df.dropna()

    #Calculate the confidence intervals for the estimates with low + high VL
    for name, group in all_estimate_df.groupby(['Group', 'Sim_Rho']):
        lower_conf = np.quantile(group['Estimated_Rho'], 0.025)
        mid_est = np.quantile(group['Estimated_Rho'], 0.5)
        upper_conf = np.quantile(group['Estimated_Rho'], 0.975)

        #append the confidence interval
        all_conf_ints.append([lower_conf, mid_est, upper_conf,
                               name[0], name[1]])

    all_conf_ints = pd.DataFrame(all_conf_ints, 
        columns=['lower_conf', 'est_rho', 'upper_conf', 'Group', 'Sim_Rho'])
    all_conf_ints.to_pickle(outDir + "all_conf_ints_" + \
                               str(NUM_BOOTSTRAPS) + "_" + str(NUM_GROUPS) \
                                + ".pkl")
    
    return


def vNE_make_comparison_dataframes(dataDir, outDir, NUM_REPS, NUM_GROUPS,
                                NUM_BOOTSTRAPS, DIST_TIME_MAX, 
                                THRESHOLD = 0.2, DSTATISTIC = False,
                                DOWNSAMPLE = True, LOCI_PER_GROUP = None,
                                PREMADE_DF = None):
    """This function randomly pairs simulations and makes a dataframe
    that shows how often we can correctly order pairs and how often
    their confidence intervals overlap. It is specifically designed for the 
    variable Ne simulations and marks the value of NE in the dataframe.
    ---------------------------------------------------------------------------
    Parameters:
    -----------
    dataDir:    str, the directory holding the data we will use for estimation.
    outDir:     str, the directory for saving the output dataframes.
    NUM_REPS:   int, the number of reps to divide into groups.
    NUM_GROUPS: int, the number of groups we will use for the analysis.
    NUM_BOOTSTRAPS: int, the number of bootstraps to use for the analysis.
    DIST_TIME_MAX:  int, the maximum distance * time to fit the curves up to.
    THRESHOLD:  float, the minimum D' value at timepoint 1 to include in the
                    analysis.
    DSTATISTIC: bool, whether or not to use the D statistic rather than D'
    DOWNSAMPLE: bool, whether or not to downsample the datasets.
    LOCI_PER_GROUP:   int, the number of segregating loci to downsample each
                        group to. Only used if DOWNSAMPLE is True.
    PREMADE_DF: str, filepath for existing d_stat dataframe to estimate from.

    Returns:
    ------------
    None, but saves the estimate and confidence interval dataframes to
         the outDir.
    """
    if not PREMADE_DF:
        #First I am going to read the data and randomly pair simulations
        all_stat_dfs = []

        #loop through each of the dataframes for the separate simulations
        for curr_data in os.listdir(dataDir):
            #only get the data directories, not hidden files
            if curr_data[0] == '.':
                continue

            #get the information for the current run
            run_info = curr_data.split('_')
            ne_val = run_info[2]
            ne_val = ne_val[2:]
            sim_rho = run_info[1]
            sim_rho = sim_rho[3:]
            rep = run_info[-1]

            if int(rep[3:]) >= NUM_REPS:
                continue

            #get the dataframe for the current run
            if DSTATISTIC:
                d_ratio_file = dataDir + curr_data + "/linkage_D/d_ratio_three_haps"
            else:
                d_ratio_file = dataDir + curr_data + "/linkage/d_ratio"

            stat_df = pd.read_pickle(d_ratio_file)

            stat_df['rep'] = int(rep[3:])
            stat_df['Ne'] = ne_val
            stat_df['Sim_Rho'] = sim_rho
            all_stat_dfs.append(stat_df)
        all_stat_dfs = pd.concat(all_stat_dfs)
    else:
        all_stat_dfs = pd.read_pickle(PREMADE_DF)

    #Filter the groups so the first d' value is greater than 0.2
    all_stat_dfs = all_stat_dfs[all_stat_dfs['d_i'] > THRESHOLD]

    #Randomly divide the reps into 10 groups
    rep_groups = np.array(range(0, NUM_REPS))
    np.random.shuffle(rep_groups)
    rep_groups = np.array_split(rep_groups, NUM_GROUPS)


    #Make a dictionary to label each group
    group_dict = {}
    for i in range(len(rep_groups)):
        for j in rep_groups[i]:
            group_dict[j] = i
            
    group_labels = [group_dict[x] for x in all_stat_dfs['rep']]
    all_stat_dfs['iter_group'] = group_labels

    print("Loaded the data Successfully")

    ########################## Estimating recombination rates ##################
    #loop through each of the distance cutoffs
    all_stat_dfs = all_stat_dfs[all_stat_dfs['Dist_X_Time'] <= DIST_TIME_MAX]
    all_stat_dfs.to_pickle(outDir + "all_stat_dfs_" + \
                               str(NUM_BOOTSTRAPS) + "_" + str(NUM_GROUPS) \
                                + ".pkl")

    all_estimate_df = []

    #loop through each rho value
    for curr_rho in all_stat_dfs['Sim_Rho'].unique():
        print(curr_rho)
        curr_rho_stat = all_stat_dfs[all_stat_dfs['Sim_Rho'] == curr_rho]

        for curr_NE in curr_rho_stat['Ne'].unique():
            curr_ne_stat = curr_rho_stat[curr_rho_stat['Ne'] == curr_NE]

            #loop through each iter group
            for curr_iteration in range(0, NUM_GROUPS):
                #get the data for the current rho and iteration
                curr_stat_df = curr_ne_stat[
                    curr_ne_stat['iter_group'] == curr_iteration]

                #downsample the group to a set number of segregating loci
                if DOWNSAMPLE:
                    curr_stat_df = plne.downsample_loci(curr_stat_df, 
                                                        LOCI_PER_GROUP)

        
                lower_fit, upper_fit, estimate_df = \
                    plne.bootstrap_rho(curr_stat_df, NUM_BOOTSTRAPS)

                estimate_df['Group'] = curr_iteration
                estimate_df['Sim_Rho'] = curr_rho
                estimate_df['Ne'] = curr_NE
                all_estimate_df.append(estimate_df)



    all_estimate_df = pd.concat(all_estimate_df, ignore_index=True)
    all_estimate_df.to_pickle(outDir + "all_estimate_df_" + \
                               str(NUM_BOOTSTRAPS) + "_" + str(NUM_GROUPS)\
                                + ".pkl")

    all_conf_ints = []

    #Drop NA values before we get the confidence interval so they don't interfere
    all_estimate_df = all_estimate_df.dropna()

    #Calculate the confidence intervals for the estimates with low + high VL
    for name, group in all_estimate_df.groupby(['Group', 'Sim_Rho', 'Ne']):
        lower_conf = np.quantile(group['Estimated_Rho'], 0.025)
        mid_est = np.quantile(group['Estimated_Rho'], 0.5)
        upper_conf = np.quantile(group['Estimated_Rho'], 0.975)

        #append the confidence interval
        all_conf_ints.append([lower_conf, mid_est, upper_conf,
                               name[0], name[1], name[2]])

    all_conf_ints = pd.DataFrame(all_conf_ints, 
        columns=['lower_conf', 'est_rho', 'upper_conf', 'Group', 'Sim_Rho', 'Ne'])
    all_conf_ints.to_pickle(outDir + "all_conf_ints_" + \
                               str(NUM_BOOTSTRAPS) + "_" + str(NUM_GROUPS) \
                                + ".pkl")
    
    return

def mig_make_comparison_dataframes(dataDir, outDir, NUM_REPS, NUM_GROUPS,
                                NUM_BOOTSTRAPS, DIST_TIME_MAX, THRESHOLD = 0.2,
                                DSTATISTIC = False, DOWNSAMPLE = True,
                                LOCI_PER_GROUP = None, PREMADE_DF = None):
    """This function randomly pairs simulations and makes a dataframe
    that shows how often we can correctly order pairs and how often
    their confidence intervals overlap. It is specifically designed for the 
    substructure simulations and marks the migration rate in the dataframe.
    ---------------------------------------------------------------------------
    Parameters:
    -----------
    dataDir:    str, the directory holding the data we will use for estimation.
    outDir:     str, the directory for saving the output dataframes.
    NUM_REPS:   int, the number of reps to divide into groups.
    NUM_GROUPS: int, the number of groups we will use for the analysis.
    NUM_BOOTSTRAPS: int, the number of bootstraps to use for the analysis.
    DIST_TIME_MAX:  int, the maximum distance * time to fit the curves up to.
    THRESHOLD:  float, the minimum D' value at timepoint 1 to include in the
                    analysis.
    DSTATISTIC: bool, whether or not to use the D statistic rather than D'
    DOWNSAMPLE: bool, whether or not to downsample the datasets.
    LOCI_PER_GROUP:   int, the number of segregating loci to downsample each
                        group to. Only used if DOWNSAMPLE is True.
    PREMADE_DF: str, filepath for existing d_stat dataframe to estimate from.
    
    Returns:
    ------------
    None, but saves the estimate and confidence interval dataframes to
         the outDir.
    """
    if not PREMADE_DF:    
        #First I am going to read the data and randomly pair simulations
        all_stat_dfs = []

        #loop through each of the dataframes for the separate simulations
        for curr_data in os.listdir(dataDir):
            #only get the data directories, not hidden files
            if curr_data[0] == '.':
                continue

            #get the information for the current run
            run_info = curr_data.split('_')
            ne_val = run_info[2]
            ne_val = ne_val[2:]
            sim_rho = run_info[1]
            sim_rho = sim_rho[3:]
            rep = run_info[-2]
            m_rate = run_info[-1]
            m_rate = m_rate[2:]

            #check if we want to include this data in the analysis
            if int(rep[3:]) >= NUM_REPS:
                continue


            #get the dataframe for the current run
            if DSTATISTIC:
                d_ratio_file = dataDir + curr_data + "/linkage_D/d_ratio_three_haps"
            else:
                d_ratio_file = dataDir + curr_data + "/linkage/d_ratio"

            stat_df = pd.read_pickle(d_ratio_file)

            stat_df['rep'] = int(rep[3:])
            stat_df['m_rate'] = m_rate
            stat_df['Sim_Rho'] = sim_rho
            all_stat_dfs.append(stat_df)
        all_stat_dfs = pd.concat(all_stat_dfs)
    else:
        all_stat_dfs = pd.read_pickle(PREMADE_DF)    

    #Filter the groups so the first d' value is greater than 0.2
    all_stat_dfs = all_stat_dfs[all_stat_dfs['d_i'] > THRESHOLD]

    #Randomly divide the reps into 10 groups
    rep_groups = np.array(range(0, NUM_REPS))
    np.random.shuffle(rep_groups)
    rep_groups = np.array_split(rep_groups, NUM_GROUPS)


    #Make a dictionary to label each group
    group_dict = {}
    for i in range(len(rep_groups)):
        for j in rep_groups[i]:
            group_dict[j] = i
            
    group_labels = [group_dict[x] for x in all_stat_dfs['rep']]
    all_stat_dfs['iter_group'] = group_labels

    print("Loaded the data Successfully")

    ########################## Estimating recombination rates ##################
    #loop through each of the distance cutoffs
    all_stat_dfs = all_stat_dfs[all_stat_dfs['Dist_X_Time'] <= DIST_TIME_MAX]
    all_stat_dfs.to_pickle(outDir + "all_stat_dfs_" + \
                               str(NUM_BOOTSTRAPS) + "_" + str(NUM_GROUPS) \
                                + ".pkl")

    all_estimate_df = []

    #loop through each rho value
    for curr_rho in all_stat_dfs['Sim_Rho'].unique():
        print(curr_rho)
        curr_rho_stat = all_stat_dfs[all_stat_dfs['Sim_Rho'] == curr_rho]

        for curr_m in curr_rho_stat['m_rate'].unique():
            curr_m_stat = curr_rho_stat[curr_rho_stat['m_rate'] == curr_m]

            #loop through each iter group
            for curr_iteration in range(0, NUM_GROUPS):
                #get the data for the current rho and iteration
                curr_stat_df = curr_m_stat[
                    curr_m_stat['iter_group'] == curr_iteration]

                #downsample the group to a set number of segregating loci
                if DOWNSAMPLE:
                    curr_stat_df = plne.downsample_loci(curr_stat_df, 
                                                        LOCI_PER_GROUP)

        
                lower_fit, upper_fit, estimate_df = \
                    plne.bootstrap_rho(curr_stat_df, NUM_BOOTSTRAPS)

                estimate_df['Group'] = curr_iteration
                estimate_df['Sim_Rho'] = curr_rho
                estimate_df['m_rate'] = curr_m
                all_estimate_df.append(estimate_df)



    all_estimate_df = pd.concat(all_estimate_df, ignore_index=True)
    all_estimate_df.to_pickle(outDir + "all_estimate_df_" + \
                               str(NUM_BOOTSTRAPS) + "_" + str(NUM_GROUPS)\
                                + ".pkl")

    all_conf_ints = []

    #Drop NA values before we get the confidence interval so they don't interfere
    all_estimate_df = all_estimate_df.dropna()

    #Calculate the confidence intervals for the estimates with low + high VL
    for name, group in all_estimate_df.groupby(['Group', 'Sim_Rho', 'm_rate']):
        lower_conf = np.quantile(group['Estimated_Rho'], 0.025)
        mid_est = np.quantile(group['Estimated_Rho'], 0.5)
        upper_conf = np.quantile(group['Estimated_Rho'], 0.975)

        #append the confidence interval
        all_conf_ints.append([lower_conf, mid_est, upper_conf,
                               name[0], name[1], name[2]])

    all_conf_ints = pd.DataFrame(all_conf_ints, 
        columns=['lower_conf', 'est_rho', 'upper_conf', 'Group', 'Sim_Rho', 'm_rate'])
    all_conf_ints.to_pickle(outDir + "all_conf_ints_" + \
                               str(NUM_BOOTSTRAPS) + "_" + str(NUM_GROUPS) \
                                + ".pkl")
    
    return

def samp_make_comparison_dataframes(dataDir, outDir, NUM_REPS, NUM_GROUPS,
                                NUM_BOOTSTRAPS, DIST_TIME_MAX, THRESHOLD = 0.2,
                                DSTATISTIC = False, DOWNSAMPLE = True,
                                LOCI_PER_GROUP = None, PREMADE_DF = None):
    """This function randomly pairs simulations and makes a dataframe
    that shows how often we can correctly order pairs and how often
    their confidence intervals overlap. It is specifically designed for the 
    sampling depth simulations and marks the dataframe with sample depth.
    ---------------------------------------------------------------------------
    Parameters:
    -----------
    dataDir:    str, the directory holding the data we will use for estimation.
    outDir:     str, the directory for saving the output dataframes.
    NUM_REPS:   int, the number of reps to divide into groups.
    NUM_GROUPS: int, the number of groups we will use for the analysis.
    NUM_BOOTSTRAPS: int, the number of bootstraps to use for the analysis.
    DIST_TIME_MAX:  int, the maximum distance * time to fit the curves up to.
    THRESHOLD:  float, the minimum D' value at timepoint 1 to include in the
                    analysis.
    DSTATISTIC: bool, whether or not to use the D statistic rather than D'
    DOWNSAMPLE: bool, whether or not to downsample the datasets.
    LOCI_PER_GROUP:   int, the number of segregating loci to downsample each
                        group to. Only used if DOWNSAMPLE is True.
    PREMADE_DF: str, filepath for existing d_stat dataframe to estimate from.
    
    Returns:
    ------------
    None, but saves the estimate and confidence interval dataframes to
         the outDir.
    """
    if not PREMADE_DF:    
        #First I am going to read the data and randomly pair simulations
        all_stat_dfs = []

        #loop through each of the dataframes for the separate simulations
        for curr_data in os.listdir(dataDir):
            #only get the data directories, not hidden files
            if curr_data[0] == '.':
                continue

            #get the information for the current run
            run_info = curr_data.split('_')
            ne_val = run_info[2]
            ne_val = ne_val[2:]
            sim_rho = run_info[1]
            sim_rho = sim_rho[3:]
            rep = run_info[-1]
            sampling_depth = run_info[3]
            sampling_depth = sampling_depth[1:]

            #check if we want to include this data in the analysis
            if int(rep[3:]) >= NUM_REPS:
                continue


            #get the dataframe for the current run
            if DSTATISTIC:
                d_ratio_file = dataDir + curr_data + "/linkage_D/d_ratio_three_haps"
            else:
                d_ratio_file = dataDir + curr_data + "/linkage/d_ratio"

            stat_df = pd.read_pickle(d_ratio_file)

            stat_df['rep'] = int(rep[3:])
            stat_df['depth'] = sampling_depth
            stat_df['Sim_Rho'] = sim_rho
            all_stat_dfs.append(stat_df)
        all_stat_dfs = pd.concat(all_stat_dfs)
    else:
        all_stat_dfs = pd.read_pickle(PREMADE_DF)    

    #Filter the groups so the first d' value is greater than 0.2
    all_stat_dfs = all_stat_dfs[all_stat_dfs['d_i'] > THRESHOLD]

    #Randomly divide the reps into 10 groups
    rep_groups = np.array(range(0, NUM_REPS))
    np.random.shuffle(rep_groups)
    rep_groups = np.array_split(rep_groups, NUM_GROUPS)


    #Make a dictionary to label each group
    group_dict = {}
    for i in range(len(rep_groups)):
        for j in rep_groups[i]:
            group_dict[j] = i
            
    group_labels = [group_dict[x] for x in all_stat_dfs['rep']]
    all_stat_dfs['iter_group'] = group_labels

    print("Loaded the data Successfully")

    ########################## Estimating recombination rates ##################
    #loop through each of the distance cutoffs
    all_stat_dfs = all_stat_dfs[all_stat_dfs['Dist_X_Time'] <= DIST_TIME_MAX]
    all_stat_dfs.to_pickle(outDir + "all_stat_dfs_" + \
                               str(NUM_BOOTSTRAPS) + "_" + str(NUM_GROUPS) \
                                + ".pkl")

    all_estimate_df = []

    #loop through each rho value
    for curr_rho in all_stat_dfs['Sim_Rho'].unique():
        print(curr_rho)
        curr_rho_stat = all_stat_dfs[all_stat_dfs['Sim_Rho'] == curr_rho]

        for curr_depth in curr_rho_stat['depth'].unique():
            curr_m_stat = curr_rho_stat[curr_rho_stat['depth'] == curr_depth]

            #loop through each iter group
            for curr_iteration in range(0, NUM_GROUPS):
                #get the data for the current rho and iteration
                curr_stat_df = curr_m_stat[
                    curr_m_stat['iter_group'] == curr_iteration]

                #downsample the group to a set number of segregating loci
                if DOWNSAMPLE:
                    curr_stat_df = plne.downsample_loci(curr_stat_df, 
                                                        LOCI_PER_GROUP)

        
                lower_fit, upper_fit, estimate_df = \
                    plne.bootstrap_rho(curr_stat_df, NUM_BOOTSTRAPS)

                estimate_df['Group'] = curr_iteration
                estimate_df['Sim_Rho'] = curr_rho
                estimate_df['depth'] = curr_depth
                all_estimate_df.append(estimate_df)



    all_estimate_df = pd.concat(all_estimate_df, ignore_index=True)
    all_estimate_df.to_pickle(outDir + "all_estimate_df_" + \
                               str(NUM_BOOTSTRAPS) + "_" + str(NUM_GROUPS)\
                                + ".pkl")

    all_conf_ints = []

    #Drop NA values before we get the confidence interval so they don't interfere
    all_estimate_df = all_estimate_df.dropna()

    #Calculate the confidence intervals for the estimates with low + high VL
    for name, group in all_estimate_df.groupby(['Group', 'Sim_Rho', 'depth']):
        lower_conf = np.quantile(group['Estimated_Rho'], 0.025)
        mid_est = np.quantile(group['Estimated_Rho'], 0.5)
        upper_conf = np.quantile(group['Estimated_Rho'], 0.975)

        #append the confidence interval
        all_conf_ints.append([lower_conf, mid_est, upper_conf,
                               name[0], name[1], name[2]])

    all_conf_ints = pd.DataFrame(all_conf_ints, 
        columns=['lower_conf', 'est_rho', 'upper_conf', 'Group', 'Sim_Rho', 'depth'])
    all_conf_ints.to_pickle(outDir + "all_conf_ints_" + \
                               str(NUM_BOOTSTRAPS) + "_" + str(NUM_GROUPS) \
                                + ".pkl")
    
    return
