#!/usr/bin/bash
#$ -pe serial 4


echo "Start - `date`"

source /net/feder/vol1/home/evromero/anaconda3/bin/activate hiv_recomb

cd /net/feder/vol1/home/evromero/2021_hiv-rec/src/revisionFigures/supp_params/

python3 rmse.py >> /net/feder/vol1/home/evromero/2021_hiv-rec/results/logs/rmse.${JOB_ID}.log 

echo "Finish - `date`"