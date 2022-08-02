import sys
import os
import numpy as np

#Process The pairs of SNPs from Zanini et al and convert them into numpy arrays

datadir = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/snpPairs'

os.system('unzip ' + datadir + "cocounts_p6_1_F1.zip")
# #if the files need to be unzipped, uncomment these lines
# fileList = os.listdir(datadir)
# for curr_file in fileList:
#     if curr_file.split('.')[-1] == 'zip':
#         os.system('unzip ' + datadir + "/" + curr_file)
#         continue

# #loop through the files and convert them to numpy arrays
# fileList = os.listdir(datadir)

# for curr_file in fileList:
#     if curr_file[0] == '.':
#         continue
#     #make sure to only get tsv files
#     file_type = curr_file.split('.')[-1]
#     if file_type == 'tsv':
#         #reshaping the arrays like the website suggested
#         curr_arr = np.loadtxt(datadir + "/" + curr_file)
#         L = int(np.sqrt(curr_arr.shape[0]/ 36))
#         curr_mat = curr_arr.reshape((6,6,L,L))

#         np.save(datadir + "/" + curr_file.split('.')[0], arr = curr_mat)

