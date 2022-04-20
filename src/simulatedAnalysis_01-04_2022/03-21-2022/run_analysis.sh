#!/usr/bin/bash
#$ -pe serial 4


echo "Start - `date`"

source /net/feder/vol1/home/evromero/anaconda3/bin/activate hiv_recomb

cd /net/feder/vol1/home/evromero/2021_hiv-rec/src/simulatedAnalysis/03-21-2022

snakemake --rerun-incomplete --cluster-config cluster.yaml --cluster  "qsub -l mfree={cluster.mfree},h_rt={cluster.h_rt}," --cores 4 -j 8 all

echo "Finish - `date`"