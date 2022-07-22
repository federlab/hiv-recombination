import requests
import sys
import os
import json


#USE WGET to download all the files with counts of pairs of SNPS

#get the reference sequence and SNPs for each patient
patientList = ['p6']

# os.system('wget --no-check-certificate https://hiv.biozentrum.unibas.ch/download/search')

#the number of timepoints for each patient
numTimepoints = {
    'p1' : 12, 'p2' : 6, 'p3' : 10, 'p4' : 8 , 'p5' : 7,
    'p6' : 7, 'p7' : 11,'p8' : 7, 'p9' : 8, 'p10' : 9,
    'p11' : 7
}
fragments = ['_F1', '_F2', '_F3', '_F4', '_F5', '_F6']

# baseURL = 'wget --no-check-certificate https://hiv.biozentrum.unibas.ch/download/cocounts_'
baseURL = "wget --no-check-certificate \"https://drive.switch.ch/index.php/s/vTuzEhZSVfVqeFA/download?path=%2F&files=cocounts_"

#loop through the patients and download all their files
for patient in patientList:
    #get the number of timepoints
    patient_times = numTimepoints[patient]
    for currTime in range(1, patient_times + 1):
        for fragment in fragments:
            print(baseURL + patient + "_" + str(currTime) + fragment + ".zip")
            os.system(baseURL + patient + "_" + str(currTime) + fragment + ".zip\"")

# baseURL = "wget --no-check-certificate https://hiv.biozentrum.unibas.ch/download/viralLoad_"

# #code to download all of the viral load data
# #loop through the patients and download all their files
# for patient in patientList:
#     os.system(baseURL + patient + ".tsv")

