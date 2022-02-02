import sys
sys.path.append('/net/gs/vol1/home/evromero/2021_hiv-rec/bin')
import os
import r2Analysis as r2
import matplotlib.pyplot as plt

vcfFile = '/net/gs/vol1/home/evromero/2021_hiv-rec/data/drosophila_test/snplist2Rdgrp.vcf'
samFile = '/net/gs/vol1/home/evromero/2021_hiv-rec/data/drosophila_test/2R_short.sam.txt'
bamFile = '/net/gs/vol1/home/evromero/2021_hiv-rec/data/drosophila_test/2R_short.bam'

os.system('samtools view -bS ' + samFile + ' > ' + bamFile)


#now that we have our VCF files, we need to start reading them in.
vcfDF = r2.get_VCF_Loci(vcfFile)

#get our pairs of loci and the corresponding sam reads
genDF = r2.get_sam_reads(vcfDF, bamFile)

#now calculate R^2 values for our locus pairs
r2List, distList, supportList = r2.r2_all_loci(genDF)

transformedSupport = np.log(map(lambda x:x+1, supportList))

#plot our results
plt.figure()
plt.scatter(distList, r2List, s = transformedSupport)
plt.xlabel("Distance Between Loci")
plt.ylabel("R^2")
plt.savefig("drosophila2R")