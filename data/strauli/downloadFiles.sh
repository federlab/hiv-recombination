#!/usr/bin/bash
#$ -pe serial 4
#$ -o logs/saveEncodings.stdout.txt
#$ -e logs/saveEncodings.stderr.txt


echo "Start - `date`"

cd /net/feder/vol1/home/evromero/2021_hiv-rec/data/strauli/

#if the data hasn't been fetched
# prefetch --option-file PRJNA543982.txt

#then move sra folder to working directory
#from ~/ncbi/public/sra

#get the reads
while read accession; do
    echo $accession
    fastq-dump --split-files sra/$accession.sra
done < $(pwd)/PRJNA543982.txt
echo "Finish - `date`"