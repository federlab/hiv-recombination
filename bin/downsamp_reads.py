import sys
#for running on cluster
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import numpy as np
import pandas as pd
import sim_downsamp

DIST_SAMPLES = 100000

#Downsampling param based on our simulation param
NUM_GENOMES = 800

#Downsampling parameters based on the Zanini et al. sequencing scheme
MIN_READ_LEN = 400
MAX_READ_LEN = 700
HALF_READ_LEN = 300

genotype_file = snakemake.input[0]

filteredGenotypes = pd.read_pickle(genotype_file)

downsampled_genotypes = sim_downsamp.downsample_genotypes(filteredGenotypes, NUM_GENOMES, DIST_SAMPLES, MIN_READ_LEN, MAX_READ_LEN, HALF_READ_LEN)

downsampled_genotypes.to_pickle(snakemake.output[0])

