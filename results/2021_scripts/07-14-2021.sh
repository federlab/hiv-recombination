#!/usr/bin/bash
#$ -pe serial 4


echo "Start - `date`"

module load samtools/1.9

cd /net/feder/vol1/home/evromero/2021_hiv-rec/results/

python3 07-14-2021.py >> logs/07-14-2021.${JOB_ID}.log 

echo "Finish - `date`"