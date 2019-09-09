#!/bin/bash
#SBATCH -p serial_requeue
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem 8000
#SBATCH -t 0-12:00:00
#SBATCH -J admixture_ibd
#SBATCH -o log/admixture_ibd.out
#SBATCH -e log/admixture_ibd.err
#SBATCH --constraint=holyib

### need these modules for odyssey
### otherwise, plink (v1.9) and admixture (v1.3) need to be in your path
#module load centos6/0.0.1-fasrc01
#module load admixture-1.3.0-fasrc01
	### this modeul is broken on odyssey
#module load plink/1.90-fasrc01
module load R
module load bcftools

./qc.sh SPECIESS /n/holylfs/LABS/informatics/ashultz/CompPopGen/SPECIES_DATASETS/SPECIESS/vcf/
