import requests
import sys
import os
import json

patientList = ['p1','p2','p3', 'p4' ,'p5','p6', 'p7' ,'p8','p9','p10','p11']
patientList = ['p1']

baseURL = 'wget --no-check-certificate https://hiv.biozentrum.unibas.ch/download/act_'

# #loop through the patients and download all their files
# for patient in patientList:
#     os.system(baseURL + patient + ".zip")

#now unzip all our files 
dataDir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/snpCalls'

fileList = os.listdir(dataDir)
for curr_file in fileList:
    print(curr_file, file = sys.stderr)
    os.system('unzip ' + dataDir + "/" + curr_file)