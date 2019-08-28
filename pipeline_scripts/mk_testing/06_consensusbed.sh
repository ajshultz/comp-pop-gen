bedtools intersect -a corCor.vcfann.bed -b _Ccornix_clean_coverage_sites_merged.bed -wa > corCor.filtered.bed
bedtools intersect -a corCor.onlyCDS.genes.bed -b corCor.filtered.bed -wa -wb > corCor.final.bed 
cat corCor.final.bed | cut -f1,2,3,4,8 > corCor.final.clean.bed
bedtools intersect -a corCor.onlyCDS.genes.bed -b _Ccornix_clean_coverage_sites_merged.bed > corCor.callable.bed
