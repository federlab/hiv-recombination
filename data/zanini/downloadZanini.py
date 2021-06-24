import requests
import sys
import os
import json


#USE WGET to download all the files with counts of pairs of SNPS

#get the reference sequence and SNPs for each patient
patientList = ['p1','p2','p3', 'p4' ,'p5','p6', 'p7' ,'p8','p9','p10','p11']
#the number of timepoints for each patient
numTimepoints = {
    'p1' : 12, 'p2' : 6, 'p3' : 10, 'p4' : 8 , 'p5' : 7,
    'p6' : 7, 'p8' : 7, 'p9' : 8, 'p10' : 9,
    'p11' : 7
}
fragments = ['_F1', '_F2', '_F3', '_F4', '_F5', '_F6']

baseURL = 'wget --no-check-certificate https://hiv.biozentrum.unibas.ch/download/cocounts_'

#loop through the patients and download all their files
for patient in patientList:
    #get the number of timepoints
    patient_times = numTimepoints[patient]
    for currTime in range(1, patient_times + 1):
        for fragment in fragments:
            os.system(baseURL + patient + "_" + str(currTime) + fragment + ".zip")