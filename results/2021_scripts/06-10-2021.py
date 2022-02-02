import sys
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
import os
import r2Analysis as r2
import numpy as np
import pandas as pd


#My goal for this is to test the r2Analysis.py file to check that it can use 
#the paired ends of the read.
vcfDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/vcfs/'
samDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/strauli/bam/'
outDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/r2_plots_moving_ave/'

sample_file = 'participant1_day1144'


#now that we have our VCF files, we need to start reading them in.
vcfDF = r2.get_VCF_Loci(vcfDir + sample_file + '.vcf')

# if there are variants
if vcfDF is not None :
    #get our pairs of loci and the corresponding sam reads
    genDF = r2.get_sam_reads(vcfDF, samDir + sample_file + '.sorted.bam')

    #now calculate R^2 values for our locus pairs
    r2List, distList, supportList = r2.r2_all_loci(genDF)