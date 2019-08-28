mv genes.gff.gz corCor.genes.gff.gz
gunzip corCor.genes.gff.gz 
grep "CDS" corCor.genes.gff > corCor.onlyCDS.gff
convert2bed -i gff < corCor.onlyCDS.gff > corCor.onlyCDS.gff.bed
cat corCor.onlyCDS.gff.bed | python corCor.genenames.py > corCor.onlyCDS.genes.bed
convert2bed -i vcf < corCor.clean.ann.vcf > corCor.clean.ann.bed
