#!/usr/bin/bash
#$ -pe serial 4


echo "Start - `date`"

source /net/feder/vol1/home/evromero/anaconda3/bin/activate hiv_recomb

cd /net/feder/vol1/home/evromero/2021_hiv-rec/src/revisionFigures/sim_fig_estimation/

python3 neutral_estimation.py >> ../../../results/logs/neutral_estimation.py.${JOB_ID}.log 

echo "Finish - `date`"