import os
import sys


dirPath = "/net/gs/vol1/home/evromero/2021_hiv-rec/data/strauli/bam/"
dirList = os.listdir(dirPath)

#loop through the files in the directory
for filename in dirList:
    if (filename.split('.'))[1] == 'bam':
        currPath = dirPath + "/"+ filename
        sortedPath = dirPath + '/' + filename.split('.')[0] + '.sorted.bam'
        os.system('samtools sort ' + currPath + ' -o ' + sortedPath )
        os.system('samtools index ' + sortedPath)