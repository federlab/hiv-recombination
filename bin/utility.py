import os
import numpy as np

def getPatientList(dir):
    """ Takes a directory of files with patients and time points. 
    Generates a dictionary where the keys are patients and the values
    are lists of the corresponding data files.
    """
    participant_files = {}
    
    #loop through all of the vcf files and group them by participant
    for file in os.listdir(dir):
        #get the participant and the date
        parDate = file.split('.')[0]
        #get the participant
        curr_participant = parDate.split('_')[0]

        #if we have seen files for this participant already add its files to entry
        if curr_participant in participant_files.keys():
            participant_files[curr_participant] += [parDate]
        else:
            participant_files[curr_participant] = [parDate]
    return participant_files