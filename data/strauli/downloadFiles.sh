#!/usr/bin/bash
#$ -pe serial 4
#$ -o logs/saveEncodings.stdout.txt
#$ -e logs/saveEncodings.stderr.txt


echo "Start - `date`"

cd /net/gs/vol1/home/evromero/2021_hiv-rec/data/strauli/

#if the data hasn't been fetched
#prefetch --option-file PRJNA543982.txt

#get the reads
while read accession; do
    echo $accession
    fastq-dump --split-files SRA_files/$accession.sra
done < $(pwd)/PRJNA543982.txt
echo "Finish - `date`"