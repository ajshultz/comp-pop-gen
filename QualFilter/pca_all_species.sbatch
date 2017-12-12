#!/bin/bash

#SBATCH -p serial_requeue
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem 24000
#SBATCH -t 0-01:00:00
#SBATCH -J angsd_qc_plot
#SBATCH -o logs/qc_plot_%j.out
#SBATCH -e logs/qc_plot_%j.err
#SBATCH --constraint=holyib

module load R/3.4.1-fasrc01
export R_LIBS_USER=$HOME/sw/R:$R_LIBS_USER

while read line
do
echo Rscript --vanilla comp-pop-gen/QualFilter/quick_pca.R ${line}
done < species_list
