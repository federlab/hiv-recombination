#!/usr/bin/bash
#$ -pe serial 4


echo "Start - `date`"

module load samtools/1.9

source /net/feder/vol1/home/evromero/anaconda3/bin/activate hiv_recomb

cd /net/feder/vol1/home/evromero/2021_hiv-rec/results/

python3 11-25-2021.py >> logs/11-25-2021.${JOB_ID}.log 

echo "Finish - `date`"