import sys
import os
from Bio.Align.Applications import MafftCommandline

parList = ["2","3","4","5","6","7","8","9","10"]


for participant in parList:
    dataDir1 = "/net/gs/vol1/home/evromero/2021_hiv-rec/data/strauli/fastq/Participant\ " + participant
    dataDir2 = "/net/gs/vol1/home/evromero/2021_hiv-rec/data/strauli/fastq/Participant " + participant
    refDir = "/net/gs/vol1/home/evromero/2021_hiv-rec/data/strauli/fastq/Reference/hiv"
    outDir = "/net/gs/vol1/home/evromero/2021_hiv-rec/data/strauli/fastq/alignments/"

    # os.system('cd ' + dataDir1)
    # os.system('ls')

    filesList = os.listdir(dataDir2)
    # print(filesList, file = sys.stderr)
    # outputFile = "participant" + participant + "alignment.fasta"
    # print(outputFile, file = sys.stderr)


    for currFile in filesList:
        stringOfPaths = ""
        fileContents = currFile.split("_")

        #if the reads are labelled 
        if len(fileContents):
            stringOfPaths += dataDir1 + "/" + currFile + " "
            
            date = fileContents[3]

            commandStr = "bwa mem -M -t 4 " + refDir + \
                        " " + stringOfPaths + "> " +  outDir + "participant" + participant + "_" + date + ".sam"

            print(commandStr, file = sys.stderr)

            os.system(commandStr)
