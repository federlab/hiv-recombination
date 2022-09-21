import os

#This script is used to compress the numpy array data. It will take in the
#directory produced by the pipeline and compress the numpy directories

def compress_numpy(pipeline_dir):
    #loop throught the subdirectories
    for curr_dir in os.listdir(pipeline_dir):
        print(curr_dir)
        #ignore hidden files
        if curr_dir[0] == '.':
            continue
        #ignore non-directories
        if not os.path.isdir(os.path.join(pipeline_dir, curr_dir)):
            continue

        #get the numpy directory that we will compress
        curr_numpy_dir = os.path.join(pipeline_dir, curr_dir, 'numpy')

        if not os.path.exists(curr_numpy_dir + ".tar.gz")):
            os.system('tar -zcvf ' + curr_numpy_dir + '.tar.gz ' + curr_numpy_dir)
