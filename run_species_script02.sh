#!/bin/bash

#Run all species through script 2

for SPECIES in Ccornix Corientalis Cbrachyrhynchos Cmonedula Cpectoralis Cfrugilegus Cdauuricus Pcollybita Ptristis Fhypoleuca Fspeculigera Pmajor Pmontanus Pdomesticus Pitaliae Phispaniolensis Tbichenovii Vchrysoptera Vcyanoptera
do
python ../comp-pop-gen/pipeline_scripts/02_download_qc.py --config ../comp-pop-gen/config_files/${SPECIES}_config.txt >> python_logs/${SPECIES}.log
done