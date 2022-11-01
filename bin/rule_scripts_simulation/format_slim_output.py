import sys
#for running on cluster
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import numpy as np
import pandas as pd
import slim_out_parse as sop
import zaniniUtil as zu

#The cutoffs for processing the segregating loci
SEG_CUTOFF = 0.01

#These parameters are used when processing the data
#10-07-2022 Elena commented out since we are going to go with the less stringent filtering
# CUTOFF = 0.03
# SUCCESS = 0.01

# This script is going to take in the output from a slim run. It will then make
# a segregating loci dataframe and a haplotype dataframe. These dataframes
# can be used as input for the recombination estimation snakemake pipeline.
outdataDir = snakemake.output[0]
print(outdataDir, file = sys.stderr)
outdataDir = outdataDir.split('/')[:-2]
outdataDir = "/".join(outdataDir)

# Get the slim output we need to parse
inputFile = snakemake.input[0]
hxb2 = "/net/feder/vol1/project/hiv_recombination/data/reference/hxb2/hxb2_env_full.fa.txt"

#Read in the hxb2 sequence we will use to get the relevant nucleotides
with open(hxb2) as f:
    hxb2_seq = f.readlines()[1]
    hxb2_seq = hxb2_seq.strip()

#Read in the slim output
with open(inputFile) as f:
    lines = f.read()
    #Separate each output section and get rid of the script text
    lines = lines.split('#OUT')[1:]

#Make dataframes to store the segregating and fixed loci output
seg_loci_slim_df = []
fixed_loci_slim_df = []

#Look at the number of reads in the data
samplenum = []

#Make a dictionary of the haplotype dictionaries keyed by timepoint
time_hap_dict = {}

#Now separately parse each section
for curr_section in lines:
    #Get the timepoint
    first_line = curr_section.split('\n', 1)[0]
    first_line = first_line.split(' ')

    #Isolate the timepoint and convert to an int
    time = int(first_line[1])
    segregating = (first_line[2] == 'SS')
    if segregating: samplenum.append(int(first_line[4]))

    #Check if there is this section includes genotype info
    contains_genomes = len(curr_section.split('Genomes:')) > 1
    
    #Separate out and process the genotype info if it exists
    if contains_genomes:
        #Separate out the genotype info
        split_results = curr_section.split('Genomes:')
        genome_data = split_results[1]
        mut_data = split_results[0]

        #Organize the genotype data
        time_hap_dict[time] = sop.parse_genome_data(genome_data)
    else: mut_data = curr_section
    
    #Put the mutation data into a nice dataframe format
    if segregating:
        seg_loci_slim_df.append(sop.parse_output_section(mut_data))
    else:
        fixed_loci_slim_df.append(sop.parse_output_section(mut_data))

#Combine the segregating and fixed loci dataframes
seg_loci_slim_df = pd.concat(seg_loci_slim_df, ignore_index=True)
fixed_loci_slim_df = pd.concat(fixed_loci_slim_df, ignore_index=True)

#Get the number of samples
samplenum = np.unique(samplenum)
if len(samplenum) > 1:
    raise ValueError('The number of sampled individuals doesn\'t\
        match across timepoints. Frequencies cannot be properly calculated.')
samplenum = samplenum[0]
    
#Initialize the dataframes we will output
all_seg_loci_df = []
all_haplotype_df = []


#Now loop through timepoints in order and make the dataframes we need to output
all_fixed_times = fixed_loci_slim_df['time'].unique()
all_seg_times = seg_loci_slim_df['time'].unique()
all_times = np.unique(np.concatenate((all_fixed_times, all_seg_times)))
all_times.sort()

for curr_time in all_times:
    #Get the fixed mutations at this timepoint
    curr_fixed_df = fixed_loci_slim_df[fixed_loci_slim_df['time'] == curr_time]

    #Loop through the mutations and update the hxb2 sequence
    for index, row in curr_fixed_df.iterrows():
        #Get the nucleotide that is mutated to
        allele = row['allele']
        #Get the position of the mutation
        mut_pos = row['position']
        #Update the hxb2 sequence
        hxb2_seq = hxb2_seq[:mut_pos] + allele + hxb2_seq[mut_pos+1:]

    #Get the segregating mutations at this timepoint
    curr_seg_df = seg_loci_slim_df[seg_loci_slim_df['time'] == curr_time]
    

    #Make the segregating loci dataframe
    print('Making segregating loci dataframe for timepoint ' + str(curr_time))
    seg_loci_df = sop.make_seg_loci_df(curr_seg_df, hxb2_seq, samplenum, SEG_CUTOFF)
    seg_loci_df['timepoint'] = curr_time
    seg_loci_df['frag_len'] = len(hxb2_seq)

    #Now we can make the haplotypes dataframe
    print('Making haplotype dataframe for timepoint ' + str(curr_time))
    relevant_genomes = time_hap_dict[curr_time]
    haplotype_df = sop.make_haplotype_df(relevant_genomes, curr_seg_df, hxb2_seq)
    haplotype_df['timepoint'] = curr_time
    

    #10-07-2022 Elena commented out since we are going to go with the less stringent filtering
    # #Lastly, We will filter the haplotype df and add it to the list
    # print('Filtering dataframes for timepoint ' + str(curr_time))
    # haplotype_df = zu.filter_genotype_df(haplotype_df, seg_loci_df, CUTOFF, SUCCESS)


    #Add our dataframes to the lists
    all_seg_loci_df.append(seg_loci_df)
    all_haplotype_df.append(haplotype_df)

#Write the dataframes to files
all_seg_loci_df = pd.concat(all_seg_loci_df, ignore_index=True)
all_seg_loci_df.to_pickle(outdataDir + '/analysis/FilteredLoci')

all_haplotype_df = pd.concat(all_haplotype_df, ignore_index=True)
all_haplotype_df.to_pickle(outdataDir + '/analysis/FilteredGenotypes')