#!/bin/bash

#SBATCH -p serial_requeue
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem 8000
#SBATCH -t 0-08:00:00
#SBATCH -J IndFilt
#SBATCH -o logs/IndFilt_%A_%a.out
#SBATCH -e logs/IndFilt_%A_%a.err
#SBATCH --constraint=holyib

#Use user-installed version of vcftools (v. 0.1.15)
module load zlib/1.2.11-fasrc01
module load tabix/0.2.6-fasrc01
source /n/home13/ashultz/sw/progs/vcftools-0.1.15/setup.sh

SPECIES=$1
SPECIES_DIR=$2
SAMPLE_FILE_DIR=$3

#Create vcf clean samples dir if it doesn't exist, and vcftools logs dir if it doesn't exist
mkdir -p ${SPECIES_DIR}/${SPECIES}/vcf_clean_samples
mkdir -p logs_vcftools

#Filter individuals due to relatedness, coverage or structure
vcftools --gzvcf ${SPECIES_DIR}/${SPECIES}/vcf/${SPECIES}_hardfilters.${SLURM_ARRAY_TASK_ID}.vcf.gz \
--remove ${SAMPLE_FILE_DIR}/${SPECIES}_samples_to_remove.txt \
--recode --recode-INFO-all \
--out ${SPECIES_DIR}/${SPECIES}/vcf_clean_samples/${SPECIES}_hardfilters_samplefiltered.${SLURM_ARRAY_TASK_ID} \
2> logs_vcftools/${SPECIES}s_remove-filtered.${SLURM_ARRAY_TASK_ID}.log

#bgzip and tabix output
bgzip ${SPECIES_DIR}/${SPECIES}/vcf_clean_samples/${SPECIES}_hardfilters_samplefiltered.${SLURM_ARRAY_TASK_ID}.recode.vcf
tabix -p vcf ${SPECIES_DIR}/${SPECIES}/vcf_clean_samples/${SPECIES}_hardfilters_samplefiltered.${SLURM_ARRAY_TASK_ID}.recode.vcf.gz


