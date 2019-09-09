#!/bin/bash
#SBATCH -p shared
#SBATCH -n 2
#SBATCH -N 1
#SBATCH -J gg_30X
#SBATCH -o ./Herato_30X_HC/logs/gg_30X_%A_%a.out
#SBATCH -e ./Herato_30X_HC/logs/gg_30X_%A_%a.err

module load jdk/1.8.0_45-fasrc01

MEM=$1

gatk --java-options "-Xmx${MEM}g -XX:ParallelGCThreads=1" GenomicsDBImport -R ./Herato_30X_HC/genome/Herato.fa --variant ./Herato_30X_HC/gvcf/SAMN08049966.30X.${SLURM_ARRAY_TASK_ID}.g.vcf.gz --variant ./Herato_30X_HC/gvcf/SAMN05224117.30X.${SLURM_ARRAY_TASK_ID}.g.vcf.gz --variant ./Herato_30X_HC/gvcf/SAMN05224208.30X.${SLURM_ARRAY_TASK_ID}.g.vcf.gz --variant ./Herato_30X_HC/gvcf/SAMN05224196.30X.${SLURM_ARRAY_TASK_ID}.g.vcf.gz --variant ./Herato_30X_HC/gvcf/SAMN05224197.30X.${SLURM_ARRAY_TASK_ID}.g.vcf.gz --variant ./Herato_30X_HC/gvcf/SAMN05224099.30X.${SLURM_ARRAY_TASK_ID}.g.vcf.gz --variant ./Herato_30X_HC/gvcf/SAMN05224098.30X.${SLURM_ARRAY_TASK_ID}.g.vcf.gz --variant ./Herato_30X_HC/gvcf/SAMN05224097.30X.${SLURM_ARRAY_TASK_ID}.g.vcf.gz --variant ./Herato_30X_HC/gvcf/SAMN05224210.30X.${SLURM_ARRAY_TASK_ID}.g.vcf.gz --variant ./Herato_30X_HC/gvcf/SAMN05224198.30X.${SLURM_ARRAY_TASK_ID}.g.vcf.gz --intervals ./Herato_30X_HC/genome/15_splits_interval_lists/Herato_${SLURM_ARRAY_TASK_ID}.interval_list --genomicsdb-workspace-path Herato/genomics_db/interval_${SLURM_ARRAY_TASK_ID}

gatk --java-options "-Xmx${MEM}g -XX:ParallelGCThreads=1" GenotypeGVCFs -R ./Herato_30X_HC/genome/Herato.fa -V gendb://./Herato_30X_HC/genomics_db/interval_${SLURM_ARRAY_TASK_ID} -O ./Herato_30X_HC/vcf/Herato.30X.${SLURM_ARRAY_TASK_ID}.vcf.gz --heterozygosity 0.025 --intervals ./Herato_30X_HC/genome/15_splits_interval_lists/Herato_${SLURM_ARRAY_TASK_ID}.interval_list

gatk --java-options "-Xmx${MEM}g -XX:ParallelGCThreads=1" VariantsToTable -V ./Herato_30X_HC/vcf/Herato.30X.${SLURM_ARRAY_TASK_ID}.vcf.gz -O ./Herato_30X_HC/stats/Herato_30X_${SLURM_ARRAY_TASK_ID}_unfilteredVCFstats.txt -F CHROM -F POS -F TYPE -F HET -F HOM-REF -F HOM-VAR -F NO-CALL -F NCALLED -F QD -F MQ -F FS -F SOR -F MQRankSum -F ReadPosRankSum -R ./Herato_30X_HC/genome/Herato.fa --intervals ./Herato_30X_HC/genome/15_splits_interval_lists/Herato_${SLURM_ARRAY_TASK_ID}.interval_list
    