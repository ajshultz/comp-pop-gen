## Modules required for pipeline scripts on cluster

module load bcftools/1.5-fasrc02
module load samtools/1.5-fasrc02
module load vcftools/0.1.14-fasrc01
module load htslib/1.5-fasrc02 
module load jdk/10.0.1-fasrc01
module load gcc/7.1.0-fasrc01
module load bedops/2.4.25-fasrc01
module load bedtools2/2.26.0-fasrc01
module load python/3.6.3-fasrc01
module load R/3.5.1-fasrc01

## if using an updated version of R, will need to reload dependencies
module load gcc/8.2.0-fasrc01 openmpi/3.1.1-fasrc01 R/3.6.1-fasrc01

