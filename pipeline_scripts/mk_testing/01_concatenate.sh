bcftools concat Ccornix_hardfilters_samplefiltered.1.recode.vcf.gz Ccornix_hardfilters_samplefiltered.2.recode.vcf.gz Ccornix_hardfilters_samplefiltered.3.recode.vcf.gz Ccornix_hardfilters_samplefiltered.4.recode.vcf.gz Ccornix_hardfilters_samplefiltered.5.recode.vcf.gz Ccornix_hardfilters_samplefiltered.6.recode.vcf.gz Ccornix_hardfilters_samplefiltered.7.recode.vcf.gz  Ccornix_hardfilters_samplefiltered.8.recode.vcf.gz Ccornix_hardfilters_samplefiltered.9.recode.vcf.gz Ccornix_hardfilters_samplefiltered.10.recode.vcf.gz -O z -o corCor_concat.vcf.gz
tabix -p vcf corCor_concat.vcf.gz
vcftools --gzvcf corCor_concat.vcf.gz --remove-filtered-all --remove-indels --min-alleles 2 --max-alleles 2 --mac 1 --max-missing 0.75 --recode --recode-INFO-all --out corCor.clean
bgzip -c corCor.clean.recode.vcf > corCor.clean.vcf.gz
tabix -p vcf corCor.clean.vcf.gz
