#!/bin/bash

#Run species through scripts 01, 02 and 03, given the species abbreviation is the first argument

SPECIES=$1

python ../comp-pop-gen/pipeline_scripts/01_download_qc.py --config ../comp-pop-gen/config_files/${SPECIES}_config.txt >> python_logs/${SPECIES}.log

python ../comp-pop-gen/pipeline_scripts/02_dedup_gather_metrics.py --config ../comp-pop-gen/config_files/${SPECIES}_config.txt >> python_logs/${SPECIES}.log

python ../comp-pop-gen/pipeline_scripts/03_haplotypecalling.py --config ../comp-pop-gen/config_files/${SPECIES}_config.txt >> python_logs/${SPECIES}.log

echo ${SPECIES} is finished