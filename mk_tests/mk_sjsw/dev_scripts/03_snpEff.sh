## Annotation of VCF with calculation of variant effects on known genes with SnpEff and cleaning with SnpSift

mkdir -p data/corCor.2.ncbi/
cp /n/holylfs/LABS/informatics/ashultz/CompPopGen/SPECIES_DATASETS/Ccornix/genome/Ccornix.fa /scratch/swuitchik/CompPopGen/snpEff/data/corCor.2.ncbi/
mv Ccornix.fa sequences.fa
gzip sequences.fa
cp /n/holylfs/LABS/informatics/ashultz/CompPopGen/SPECIES_DATASETS/Ccornix/genome/GCF_000738735.2_ASM73873v2_genomic.gff.gz /scratch/swuitchik/CompPopGen/snpEff/data/corCor.2.ncbi/
mv GCF_000738735.2_ASM73873v2_genomic.gff.gz genes.gff.gz

## add the following into the snpEff.config file under the Databases & Genomes section: 

# Hooded crow genome, NCBI version 2
corCor.2.ncbi.genome : Corvus_cornix_cornix

## 

java -jar snpEff.jar build -gff3 -v corCor.2.ncbi #builds database
java -jar snpEff.jar corCor.2.ncbi data/corCor.clean.vcf.gz > corCor.ann.vcf #annotates VCF

