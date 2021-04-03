import pandas as pd
import os
import sys

#path to the directory with the Fastq files we want to sort and move
dirPath = "fastq/"

#path to the csv file with information on each run
infoFile = "SraRunInfo.csv"

#open the csv with all of the data
runInfo = pd.read_csv(infoFile)

print(runInfo, file=sys.stderr )

#iterate through the fastq files in the directory
for filename in os.listdir(dirPath):
    #just get the first version (I accidentally got duplicates of the files)
    version = filename.split('_')
    accession = version[0]
    version = version[1].split('.')
    version = int(version[0])
    if version == 1:
        print(filename, file = sys.stderr)

        #now we need to determine if the file is env sequences
        currRow = runInfo.loc[runInfo['Run'] == accession]
        #check if the file is viral RNA
        isRNA = (currRow['LibrarySource'] == 'VIRAL RNA')
        if isRNA:
            patient = currRow['SampleName']
            currLibrary = currRow['LibraryName']
            print(currLibrary, file = sys.stderr)
        
            #make sure to make the patient directory
            if not os.path.exists(patient):
                os.makedir(patient)
                
            #rename and move the file 
            newName = filename.split('.')
            newName = newName[0] + currLibrary + ".fastq"
            print(newName, file=sys.stderr)
            os.rename(dirPath + filename, dirPath + patient + "/" + newName)
        #move the antibody sequences to a different directory
        else:
            os.rename(dirPath + filename, dirPath + "transcriptome/" + newName)
        
