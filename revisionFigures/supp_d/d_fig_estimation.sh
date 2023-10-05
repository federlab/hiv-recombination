#!/usr/bin/bash
#$ -pe serial 4


echo "Start - `date`"

source /net/feder/vol1/home/evromero/anaconda3/bin/activate hiv_recomb

cd /net/feder/vol1/home/evromero/2021_hiv-rec/src/revisionFigures/supp_d/

python3 d_fig_estimation.py >> ../../../results/logs/d_fig_estimation.py.${JOB_ID}.log 

echo "Finish - `date`"