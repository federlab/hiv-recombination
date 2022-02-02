import sys
sys.path.append('/net/gs/vol1/home/evromero/2021_hiv-rec/bin')
import os
import makeVCFs as vcf

#Today I want to get vcf files to check the variants I am dealing with

dataDir = '/net/gs/vol1/home/evromero/2021_hiv-rec/data/strauli/bam/'

outDir = '/net/gs/vol1/home/evromero/2021_hiv-rec/results/vcfs'

reference = '/net/gs/vol1/home/evromero/2021_hiv-rec/data/strauli/fastq/Reference/HXB2Ref.fasta'

if not os.path.exists(outDir):
    os.mkdir(outDir)

vcf.createVCF(dataDir, reference, outDir + "/")

# bcfDir = "/net/gs/vol1/home/evromero/2021_hiv-rec/results/bcfs"

# vcf.bcfToVCF(bcfDir, outDir)