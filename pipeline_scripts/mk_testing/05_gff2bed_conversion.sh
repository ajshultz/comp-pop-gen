## Extraction of CDS regions & conversion to BED file

mv sequences.fa corCor.sequences.fa
cat corCor.sequences.fa | bcftools consensus corCor.clean.ann.vcf.gz -o vcfout.fasta
grep "CDS" corCor.genes.gff > corCor.onlyCDS.gff
convert2bed -i gff < corCor.onlyCDS.gff > corCor.onlyCDS.gff.bed
