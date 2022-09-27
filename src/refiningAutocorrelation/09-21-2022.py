import pandas as pd

#Today I am just looking at the formats of the data for the segregating loci
#and genotype DFs. I'm hoping just to put the SLIM data directly in this format.

input_dir = "/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini_snakemake/p1_F1/analysis/"

#read in the segregating loci data
#times are in generations
seg_loc_df = pd.read_pickle(input_dir + "FilteredLoci")
print(seg_loc_df.head())
print(seg_loc_df['timepoint'].unique())

genotype_df = pd.read_pickle(input_dir + "FilteredGenotypes")
print(genotype_df.head())
