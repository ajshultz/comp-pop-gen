#!/bin/bash

#SBATCH -p shared
#SBATCH -n 2
#SBATCH -N 1
#SBATCH --mem 16000
#SBATCH -t 0-8:00
#SBATCH -o ./logs/dedup_%j.out
#SBATCH -e ./logs/dedup_%j.err
#SBATCH --constraint=holyib

SAMPLE=$1

module load java/1.8.0_45-fasrc01

gatk --java-options "-Xmx8g -XX:ParallelGCThreads=1" MarkDuplicatesGATK \
-I ${SAMPLE}.sorted.bam \
-O ${SAMPLE}.dedup.gatk.bam \
--METRICS_FILE ${SAMPLE}.dedup.gatk.metrics.txt \
REMOVE_DUPLICATES false