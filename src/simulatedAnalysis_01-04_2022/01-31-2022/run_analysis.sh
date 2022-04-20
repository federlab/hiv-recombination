#!/usr/bin/bash
#$ -pe serial 4


echo "Start - `date`"

source /net/feder/vol1/home/evromero/anaconda3/bin/activate hiv_recomb

cd /net/feder/vol1/home/evromero/2021_hiv-rec/results/simulatedAnalysis/01-31-2022

snakemake --cluster "qsub -l mfree=100G,h_rt=120:0:0," --cores 4 -j 4 all

echo "Finish - `date`"