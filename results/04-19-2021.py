import sys
sys.path.append('/net/gs/vol1/home/evromero/2021_hiv-rec/bin')
import os
import r2Analysis as r2
import matplotlib.pyplot as plt
import numpy as np

vcfDir = '/net/gs/vol1/home/evromero/2021_hiv-rec/results/vcfs/'
samDir = '/net/gs/vol1/home/evromero/2021_hiv-rec/data/strauli/bam/'
outDir = '/net/gs/vol1/home/evromero/2021_hiv-rec/results/r2_plots/'

for file in os.listdir(vcfDir):
    print(file, sys.stderr)
    parDate = file.split('.')[0]

    #now that we have our VCF files, we need to start reading them in.
    vcfDF = r2.get_VCF_Loci(vcfDir + parDate + '.vcf')

    # if there are variants
    if vcfDF is not None :
        #get our pairs of loci and the corresponding sam reads
        genDF = r2.get_sam_reads(vcfDF, samDir + parDate + '.sorted.bam')

        #now calculate R^2 values for our locus pairs
        r2List, distList, supportList = r2.r2_all_loci(genDF)

        #plot our results
        plt.figure()
        plt.scatter(distList, r2List, s = supportList, alpha = 0.5)
        l1 = plt.scatter([],[], s=10, edgecolors='none')
        l2 = plt.scatter([],[], s=100, edgecolors='none')
        l3 = plt.scatter([],[], s=1000, edgecolors='none')
        l4 = plt.scatter([],[], s=5000, edgecolors='none')
        labels = ["10", "100", "1000", "5000"]
        leg = plt.legend([l1, l2, l3, l4], labels, ncol=4, frameon=True, fontsize=12,
        loc="upper center", bbox_to_anchor=(0.5, 1.5),
        handlelength=2, borderpad = 1.8,
        handletextpad=1, scatterpoints = 1)
        plt.xlabel("Distance Between Loci")
        plt.ylabel("R^2")
        plt.tight_layout()
        plt.savefig(outDir + parDate)
        plt.close()


