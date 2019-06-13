#!/bin/bash

#This script will create filtered VCFs given a species abbreviation (argument 1), path to existing SPECIES_DATASETS directory (argument 2), path to file with samples to remove (argument 3). Make sure that the array is equal to the number of VCF files for the species of interest

#Check if logs directory exists, if not, create one
mkdir -p logs

#Set common variables
SPECIES_DIR=/n/holylfs/LABS/informatics/ashultz/CompPopGen/SPECIES_DATASETS
SAMPLE_FILE_DIR=/n/holylfs/LABS/informatics/ashultz/CompPopGen/comp-pop-gen/vcf_qc/samples_to_remove

#Submit sample removal jobs 

sbatch --array=1-10 /n/holylfs/LABS/informatics/ashultz/CompPopGen/comp-pop-gen/vcf_qc/vcf_sample_filter.sbatch Acunicularia ${SPECIES_DIR} ${SAMPLE_FILE_DIR} 

sbatch --array=1-10 /n/holylfs/LABS/informatics/ashultz/CompPopGen/comp-pop-gen/vcf_qc/vcf_sample_filter.sbatch Ccornix ${SPECIES_DIR} ${SAMPLE_FILE_DIR}

sbatch --array=1-35 /n/holylfs/LABS/informatics/ashultz/CompPopGen/comp-pop-gen/vcf_qc/vcf_sample_filter.sbatch Falbicollis ${SPECIES_DIR} ${SAMPLE_FILE_DIR}

sbatch --array=1-32 /n/holylfs/LABS/informatics/ashultz/CompPopGen/comp-pop-gen/vcf_qc/vcf_sample_filter.sbatch Pdomesticus ${SPECIES_DIR} ${SAMPLE_FILE_DIR}

sbatch --array=1-34 /n/holylfs/LABS/informatics/ashultz/CompPopGen/comp-pop-gen/vcf_qc/vcf_sample_filter.sbatch Pmajor ${SPECIES_DIR} ${SAMPLE_FILE_DIR}

sbatch --array=1-10 /n/holylfs/LABS/informatics/ashultz/CompPopGen/comp-pop-gen/vcf_qc/vcf_sample_filter.sbatch Ptrochilus ${SPECIES_DIR} ${SAMPLE_FILE_DIR}