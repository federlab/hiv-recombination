import os
import sys

dirPath = "/net/gs/vol1/home/evromero/2021_hiv-rec/data/strauli/fastq/alignments"
dirList = os.listdir(dirPath)
outPath = "/net/gs/vol1/home/evromero/2021_hiv-rec/data/strauli/bam/"


#loop through the files in the directory
for filename in dirList:
    outFile = outPath + (filename.split('.'))[0] + '.bam'
    os.system('samtools view -bS ' + dirPath + "/"+ filename + ' > ' + outFile)