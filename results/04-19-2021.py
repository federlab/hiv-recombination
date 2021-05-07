import sys
sys.path.append('/net/gs/vol1/home/evromero/2021_hiv-rec/bin')
import os
import r2Analysis as r2
import matplotlib.pyplot as plt

vcfDir = '/net/gs/vol1/home/evromero/2021_hiv-rec/results/vcfs/'
samDir = '/net/gs/vol1/home/evromero/2021_hiv-rec/data/strauli/bam/'

#now that we have our VCF files, we need to start reading them in.
vcfDF = r2.get_VCF_Loci(vcfDir + 'participant1_day4903.vcf')

#get our pairs of loci and the corresponding sam reads
genDF = r2.get_sam_reads(vcfDF, samDir + 'participant1_day4903.sorted.bam')

#now calculate R^2 values for our locus pairs
r2List, distList = r2.r2_all_loci(genDF)

#plot our results
plt.figure()
plt.scatter(distList, r2List)
plt.xlabel("Distance Between Loci")
plt.ylabel("R^2")
plt.savefig("participant1_day4903-filteredFreqs")

