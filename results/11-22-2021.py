import sys
#for cluster run
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import neher
import zaniniUtil as zu

#directories for cluster run
dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/analysis/'
viralLoadDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/viralLoads/'
outData = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/analysis/neher_downsampled/'

#Today I am going to conduct the neher and leitner analysis
#But I am going to downsample multinomially at each timepoint
#I need to do this, save the data, then match up the timepoints with viral loads.
#Note: we are not filtering the haplotypes for this

#Next we need to set some things for the run
fragment_list = ['F1','F2', 'F3', 'F4', 'F5','F6']
par_list = ['p1', 'p2', 'p3', 'p5', 'p6', 'p8', 'p9', 'p10', 'p11']
# par_list =['p1', 'p2']
# fragment_list = ['F1']
available_files_hap = os.listdir(dataDir + "haplotype_dfs/")
available_files_seg = os.listdir(dataDir + "segregating_Loci/")


NUM_VL_GROUPS = 5
BINWIDTH = 100
MIN_BIN = 0
MAX_BIN = 700
CUTOFF = 0.03
#The frequency a haplotype must have to trigger a test success
SUCCESSFILTER = 0
NUM_DOWNSAMPS = 100
RUNNAME = str(CUTOFF) +  "filtering_success_"  + str(SUCCESSFILTER) +"_binnedVL.png"

# #We can start by getting the viral load data
# viralLoadData = zu.make_viral_load_df(viralLoadDir)
# #We can also load our dataframe for converting days to timepoints
# dayToTime = pd.read_csv(dataDir + "DaysToTimepoints.csv")
# labeledLoads = []

# #We need to add a timepoint column to our viral load data
# for curr_par in par_list:
#     curr_data = viralLoadData[viralLoadData['Participant'] == curr_par]
#     curr_days = dayToTime[dayToTime['Participant'] == curr_par]
#     labels = []
#     #get the timepoint for the day
#     for index, row in curr_data.iterrows():
#         date = int(row['Days from infection'])
#         currLabel = curr_days[curr_days[' Day'] == date]
#         # print("****************", file = sys.stderr)
#         # print(curr_par, file = sys.stderr)
#         # print(currLabel, file = sys.stderr)
#         # print(date, file = sys.stderr)
#         labels.append(currLabel[' Timepoint'].tolist()[0])
#     labeled_current = curr_data.copy()
#     labeled_current['Timepoint'] = labels
#     labeledLoads.append(labeled_current)
# labeledLoads = pd.concat(labeledLoads)
# viralLoadData = labeledLoads

print("Added Timepoint Column to Viral Load Data", file = sys.stderr)

#make a list to store our haplotype dataframes in
hap_dfs = []
#make a dataframe to store our data to plot
all_frequencies_patients = []

#loop through the patients and get their results from the neher analysis
for curr_par in par_list:
    for curr_fragment in fragment_list:
        par_frag = curr_par +  "_" + curr_fragment
        haplotype_file = "haplotypes_" + par_frag + ".pkl"
        loci_file = "segregating_" + par_frag + ".pkl"
        
        #check if there are  files for this
        if haplotype_file not in available_files_hap or loci_file not in available_files_seg:
            continue

        #First we need to load in the dataframes for each patient + fragment.
        haplotype_df = pd.read_pickle(dataDir + "haplotype_dfs/" + haplotype_file)
        segregating_Loci = pd.read_pickle(dataDir + "segregating_Loci/" + loci_file)

        #now we can perform neher analysis on this dataframe
        #first we will filter out genotypes with alleles at less than 3% frequency
        haplotype_df = zu.filter_genotype_df(haplotype_df, segregating_Loci, allele_cutoff = CUTOFF, hap_cutoff = SUCCESSFILTER)
        if haplotype_df.empty:
            continue


        print("Finished collecting the Haplotype Dataframes", file = sys.stderr)

        downsampled_haps = []

        print("Finished one resampling")
        #Now we want to make sure all our genotypes have frequency > hap_cutoff
        for name, group in haplotype_df.groupby(['Locus_1', 'Locus_2', 'timepoint']):
            group.reset_index(inplace = True)
            #get the number of haplotypes
            num_haps = len(group['2_Haplotype'].unique())
            #make dictionaries so we can easily associate rows with haplotypes and frequencies
            group_dict = group.to_dict()
            freq_dict = group_dict['hap_freq']
            hap_dict = group_dict['2_Haplotype']

            #make the probability vector for our outcomes
            hap_freqs = [0 for x in range(num_haps)]
            for i in range(num_haps):
                hap_freqs[i] = freq_dict[i]

            #Downsample from a multinomial
            multinom_draw = np.random.multinomial(n= 20, pvals = hap_freqs ,size = NUM_DOWNSAMPS)

            #Now, make new dataframe entries for each of the samples
            for i in range(num_haps):
                #I think we can just new column values for how many times haplotypes were sampled
                haplotype_row = group.iloc[i]
                haplotype_row = haplotype_row.copy()
                for sample in range(NUM_DOWNSAMPS):
                    haplotype_row['Sample_' + str(sample)] = multinom_draw[sample, i]
                downsampled_haps.append(haplotype_row) 

        downsampled_haps= pd.DataFrame(downsampled_haps, columns = downsampled_haps[0].index)
        all_recombs = []
        all_muts = []
        print(downsampled_haps, file = sys.stderr)

        print("Downsampled Successfully", file = sys.stderr)

        #now for each sample, calculate the 
        for sample in range(NUM_DOWNSAMPS):
            print("entered loop", file = sys.stderr)
            print(downsampled_haps, file = sys.stderr)
            curr_sample_haps = downsampled_haps[downsampled_haps['Sample_' + str(sample)].gt(0)]
            print(curr_sample_haps, file = sys.stderr)
            recombination_df, mutation_df = neher.run_analysis(curr_sample_haps, verbose = False, success_filt = SUCCESSFILTER)

            print("Ran Neher Tests", file = sys.stderr)

        # #Now get the average viral load between the two timepoints at which the tests were conducted
        # #There is probably a clever way to do this in pandas, but I can't think of it
        # average_vls = []
        # participant_vls = viralLoadData[viralLoadData['Participant'] == curr_par]
        # #I am just going to loop through the dataframe and perform each calculation one by one
        # for index, row in mutation_df.iterrows():
        #     time_bef = int(row['Last_Timepoint'])
        #     success_time = int(row['Curr_Timepoint'])

        #     #get the viral load at the timepoint that triggered the test
        #     vl_bef = participant_vls[participant_vls['Timepoint'] == time_bef]
        #     if vl_bef.empty:
        #         average_vls.append(float('inf'))
        #     else:
        #          vl_bef = vl_bef['Viral load [virions/ml]'].tolist()[0]
            
        #     #get the viral load at the success timepoint
        #     vl_success = participant_vls[participant_vls['Timepoint'] == success_time]
        #     if vl_success.empty:
        #         average_vls.append(float('inf'))
        #     else: 
        #         vl_success = vl_success['Viral load [virions/ml]'].tolist()[0]
        #         average_vls.append(sum([vl_success, vl_bef])/2)
        # mutation_df['Average_vl'] = average_vls

        # #now repeat the same thing for our recombination df to get the average viral loads
        # average_vls = []
        # for index, row in recombination_df.iterrows():
        #     time_bef = int(row['Last_Timepoint'])
        #     success_time = int(row['Curr_Timepoint'])

        #     #get the viral load at the timepoint that triggered the test
        #     vl_bef = participant_vls[participant_vls['Timepoint'] == time_bef]
        #     if vl_bef.empty:
        #         average_vls.append(float('inf'))
        #     else: vl_bef = vl_bef['Viral load [virions/ml]'].tolist()[0]

        #     #get the viral load at the success timepoint
        #     vl_success = participant_vls[participant_vls['Timepoint'] == success_time]
        #     if vl_success.empty:
        #         average_vls.append(float('inf'))
        #     else:
        #         vl_success = vl_success['Viral load [virions/ml]'].tolist()[0]
        #         average_vls.append(sum([vl_success, vl_bef])/2)
        # recombination_df['Average_vl'] = average_vls
            recombination_df['Sample'] = sample
            mutation_df['Sample'] = sample
            all_recombs.append(recombination_df)
            all_muts.append(mutation_df)
            # write our dataframes to files
            recombination_df.to_csv(outData + "Recombination_"+ curr_par + "_" + curr_fragment + "_" + str(sample) +".csv")
            mutation_df.to_csv(outData + "Mutation_"+ curr_par + "_" + curr_fragment + "_" + str(sample) + ".csv")



