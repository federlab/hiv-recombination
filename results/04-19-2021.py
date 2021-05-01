import sys
sys.path.append('/net/gs/vol1/home/evromero/2021_hiv-rec/bin')
import os
import r2Analysis as r2

vcfDir = '/net/gs/vol1/home/evromero/2021_hiv-rec/results/vcfs/'
samDir = '/net/gs/vol1/home/evromero/2021_hiv-rec/data/strauli/bam/'

#now that we have our VCF files, we need to start reading them in.
vcfDF = r2.get_VCF_Loci(vcfDir + 'participant1_day1144.vcf')

#get our pairs of loci and the corresponding sam reads
r2.get_sam_reads(vcfDF, samDir + 'participant1_day1144.sorted.bam')

