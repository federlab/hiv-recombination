import sys
import os
from Bio.Align.Applications import MafftCommandline

parList = ["1","2","3","4","5","6","7","8","9","10"]

for participant in parList:
    dataDir1 = "/net/feder/vol1/home/evromero/2021_hiv-rec/data/strauli/fastq2/Participant\ " + participant
    dataDir2 = "/net/feder/vol1/home/evromero/2021_hiv-rec/data/strauli/fastq2/Participant " + participant
    refDir = "/net/feder/vol1/home/evromero/2021_hiv-rec/data/strauli/fastq2/Reference/hiv"
    outDir = "/net/feder/vol1/home/evromero/2021_hiv-rec/data/strauli/fastq2/alignments/"

    # os.system('cd ' + dataDir1)
    # os.system('ls')

    filesList = os.listdir(dataDir2)
    #only get files
    # outputFile = "participant" + participant + "alignment.fasta"
    # print(outputFile, file = sys.stderr)


    for currFile in filesList:
        stringOfPaths = ""
        fileContents = currFile.split("_")
        print(currFile, file = sys.stderr)

        #if the reads are labelled 
        if fileContents[-1] == "labeled.fastq":
            stringOfPaths += dataDir1 + "/" + currFile + " "
            
            date = fileContents[3]
            end = fileContents[1]

            commandStr = "bwa mem -M -t 4 " + refDir + \
                        " " + stringOfPaths + "> " +  outDir + "participant" + participant + "_" + \
                             date + "_" + end +  ".sam"

            print(commandStr, file = sys.stderr)

            os.system(commandStr)
