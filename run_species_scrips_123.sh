#!/bin/bash

#Run all species through scripts 1, 2 and 3

for SPECIES in Acunicularia Cjaponica Fheterclitus
do
python ../comp-pop-gen/pipeline_scripts/01_download_qc.py --config ../comp-pop-gen/config_files/${SPECIES}_config.txt >> python_logs/${SPECIES}.log

python ../comp-pop-gen/pipeline_scripts/02_dedup_gather_metrics.py --config ../comp-pop-gen/config_files/${SPECIES}_config.txt >> python_logs/${SPECIES}.log

python ../comp-pop-gen/pipeline_scripts/03_haplotypecalling.py --config ../comp-pop-gen/config_files/${SPECIES}_config.txt >> python_logs/${SPECIES}.log
done