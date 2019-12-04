# cp genes.gff.gz and corCor.ann.vcf to a new working directory

mv genes.gff.gz corCor.genes.gff.gz
gunzip corCor.genes.gff.gz 
grep "CDS" corCor.genes.gff > corCor.onlyCDS.gff
convert2bed -i gff < corCor.onlyCDS.gff > corCor.onlyCDS.gff.bed
cat corCor.onlyCDS.gff.bed | python corCor.genenames.py > corCor.onlyCDS.genes.bed

bedtools intersect -a corCor.vcfann.bed -b _Ccornix_clean_coverage_sites_merged.bed -wa > corCor.filtered.bed
bedtools intersect -a corCor.onlyCDS.genes.bed -b corCor.filtered.bed -wa -wb > corCor.final.bed 
cat corCor.final.bed | cut -f1,2,3,4,8 > corCor.final.clean.bed
bedtools intersect -a corCor.onlyCDS.genes.bed -b _Ccornix_clean_coverage_sites_merged.bed > corCor.callable.bed
