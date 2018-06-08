#!/bin/bash
#SBATCH -p shared
#SBATCH -n 2
#SBATCH -N 1
#SBATCH -J gg_15X
#SBATCH -o ./Mmulatta_15X_HC/logs/gg_15X_%A_%a.out
#SBATCH -e ./Mmulatta_15X_HC/logs/gg_15X_%A_%a.err

module load jdk/1.8.0_45-fasrc01

MEM=$1

gatk --java-options "-Xmx${MEM}g -XX:ParallelGCThreads=1" GenomicsDBImport -R ./Mmulatta_15X_HC/genome/Mmulatta.fa --variant ./Mmulatta_15X_HC/gvcf/SAMN03083651.15X.${SLURM_ARRAY_TASK_ID}.g.vcf.gz --variant ./Mmulatta_15X_HC/gvcf/SAMN03264600.15X.${SLURM_ARRAY_TASK_ID}.g.vcf.gz --variant ./Mmulatta_15X_HC/gvcf/SAMN03264598.15X.${SLURM_ARRAY_TASK_ID}.g.vcf.gz --variant ./Mmulatta_15X_HC/gvcf/SAMN03264597.15X.${SLURM_ARRAY_TASK_ID}.g.vcf.gz --variant ./Mmulatta_15X_HC/gvcf/SAMN03264613.15X.${SLURM_ARRAY_TASK_ID}.g.vcf.gz --variant ./Mmulatta_15X_HC/gvcf/SAMN03264638.15X.${SLURM_ARRAY_TASK_ID}.g.vcf.gz --variant ./Mmulatta_15X_HC/gvcf/SAMN03264639.15X.${SLURM_ARRAY_TASK_ID}.g.vcf.gz --variant ./Mmulatta_15X_HC/gvcf/SAMN03264641.15X.${SLURM_ARRAY_TASK_ID}.g.vcf.gz --variant ./Mmulatta_15X_HC/gvcf/SAMN03264637.15X.${SLURM_ARRAY_TASK_ID}.g.vcf.gz --variant ./Mmulatta_15X_HC/gvcf/SAMN03264635.15X.${SLURM_ARRAY_TASK_ID}.g.vcf.gz -O ./Mmulatta_15X_HC/gvcf/Mmulatta.15X.${SLURM_ARRAY_TASK_ID}.g.vcf.gz --intervals ./Mmulatta_15X_HC/genome/CHROMOSOME_splits_interval_lists/Mmulatta_${SLURM_ARRAY_TASK_ID}.interval_list --genomicsdb-workspace-path ./Mmulatta_15X_HC/genomics_db/interval_${SLURM_ARRAY_TASK_ID}

gatk --java-options "-Xmx${MEM}g -XX:ParallelGCThreads=1" GenotypeGVCFs -R ./Mmulatta_15X_HC/genome/Mmulatta.fa -V gendb://./Mmulatta_15X_HC/genomics_db/interval_${SLURM_ARRAY_TASK_ID} -O ./Mmulatta_15X_HC/vcf/Mmulatta.15X.${SLURM_ARRAY_TASK_ID}.vcf.gz --heterozygosity 0.0024 --intervals ./Mmulatta_15X_HC/genome/CHROMOSOME_splits_interval_lists/Mmulatta_${SLURM_ARRAY_TASK_ID}.interval_list

gatk --java-options "-Xmx${MEM}g -XX:ParallelGCThreads=1" VariantsToTable -V ./Mmulatta_15X_HC/vcf/Mmulatta.15X.${SLURM_ARRAY_TASK_ID}.vcf.gz -O ./Mmulatta_15X_HC/stats/Mmulatta_15X_${SLURM_ARRAY_TASK_ID}_unfilteredVCFstats.txt -F CHROM -F POS -F TYPE -F HET -F HOM-REF -F HOM-VAR -F NO-CALL -F NCALLED -F QD -F MQ -F FS -F SOR -F MQRankSum -F ReadPosRankSum -R ./Mmulatta_15X_HC/genome/Mmulatta.fa --intervals ./Mmulatta_15X_HC/genome/CHROMOSOME_splits_interval_lists/Mmulatta_${SLURM_ARRAY_TASK_ID}.interval_list
    
