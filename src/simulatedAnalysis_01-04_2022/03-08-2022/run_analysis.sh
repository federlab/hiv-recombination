#!/usr/bin/bash
#$ -pe serial 4


echo "Start - `date`"

source /net/feder/vol1/home/evromero/anaconda3/bin/activate hiv_recomb

cd /net/feder/vol1/home/evromero/2021_hiv-rec/src/simulatedAnalysis/03-08-2022

snakemake -R calc_linkage --latency-wait 10 --wait-for-files --cluster  "qsub -l mfree=8G,h_rt=120:0:0," --cores 4 -j 10 

echo "Finish - `date`"