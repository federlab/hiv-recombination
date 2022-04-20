#!/usr/bin/bash
#$ -pe serial 4


echo "Start - `date`"

source /net/feder/vol1/home/evromero/anaconda3/bin/activate hiv_recomb

cd /net/feder/vol1/home/evromero/2021_hiv-rec/src/simulatedAnalysis

python3 03-14-2022.py >> /net/feder/vol1/home/evromero/2021_hiv-rec/results/logs/03-14-2022.${JOB_ID}.log 

echo "Finish - `date`"