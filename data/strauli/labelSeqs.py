import pandas as pd
import os
import sys
from Bio import SeqIO

parList = ['3', '4', '5', '6', '7', '8', '9','10' ]
for currPar in parList:
    #path to the directory with the Fastq files we want to work with
    dirPath = "/net/gs/vol1/home/evromero/2021_hiv-rec/data/strauli/fastq/Participant "+ currPar

    dirList = os.listdir(dirPath)
    #loop through the files in the directory
    for filename in dirList:
        currDate = (filename.split('_'))[3]
        print(currDate, file = sys.stderr)

        #create a name for our labeled file and open it
        output_file_name = dirPath + '/' + filename.split('.')[0] + "_labeled.fastq"
        with open(output_file_name, 'w+') as corrected:
            for record in SeqIO.parse(open(dirPath + '/' +filename),'fastq'):
                record.name = record.id + "_" + currDate
                record.description = record.id + "_" + currDate
                record.id = record.id + "_" + currDate
                
                print(record.name, file = sys.stderr)
                SeqIO.write(record, corrected, 'fastq')
