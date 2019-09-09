#!/bin/bash

### parameters to be passed via the command line are
SPECIES=$1    	### species id, not important beyond linking all output files under that id
VCF=$2		### full path to the directory containing vcf files

### check number of individuals first
n_ind=`gzip -dc ${VCF}*.vcf.gz | grep -e '#CHROM' | awk '{ print NF; exit; }'`

if [[ $n_ind < 11 ]]; then
        echo "too few individuals, found: $n_ind - 10" 
        exit
fi

### import file and label variants
### want to go through and do some LD pruning (be aggressive)
### then print out inbreeding stats  
plink --vcf <( bcftools concat $VCF/*hardfilter*.vcf.gz | grep -v '\*' ) --make-bed --out $SPECIES --set-missing-var-ids @:#[$SPECIES]\$1,\$2 --allow-extra-chr
plink --bfile ${SPECIES} --indep-pairwise 500 10 0.1 --out ${SPECIES} --allow-extra-chr
plink --bfile ${SPECIES} --make-bed --extract ${SPECIES}.prune.in --out ${SPECIES}.ld_pruned --allow-extra-chr
plink --bfile ${SPECIES}.ld_pruned --ibc --out ${SPECIES} --allow-extra-chr

### make plot of inbreeding coefficients
echo "x <- read.table(\"${SPECIES}.ibc\",header=T)" > ${SPECIES}.plot.r
echo "pdf(\"${SPECIES}.ibc.pdf\",height=5,width=5)" >> ${SPECIES}.plot.r
echo "plot(x\$Fhat1,x\$Fhat2,xlab=\"Fhat1\",ylab=\"Fhat2\")" >> ${SPECIES}.plot.r
echo "dev.off()" >> ${SPECIES}.plot.r
Rscript ${SPECIES}.plot.r

### run pca
plink -bfile ${SPECIES}.ld_pruned --pca --out ${SPECIES} --allow-extra-chr

### plot it 
echo "pdf(\"${SPECIES}.pca.pdf\",height=8,width=5)" > ${SPECIES}.pca.plot
echo "par(mfrow=c(2,1),mar=c(4,4,2,2))" >> ${SPECIES}.pca.plot
echo "d <- read.table(\"${SPECIES}.eigenval\")" >> ${SPECIES}.pca.plot
echo "plot(c(seq(1,length(d\$V1),by=1)),d\$V1/sum(d\$V1)*100,xlab=\"PC\",ylab=\"Percent Variance Explained\")" >> ${SPECIES}.pca.plot
echo "d <- read.table(\"${SPECIES}.eigenvec\")" >> ${SPECIES}.pca.plot
echo "plot(d\$V3,d\$V4,cex=0.5,xlab=\"PC 1\",ylab = \"PC 2\")" >> ${SPECIES}.pca.plot
echo "dev.off()" >> ${SPECIES}.pca.plot
Rscript ${SPECIES}.pca.plot

### pairwise IBD analysis
plink -bfile ${SPECIES}.ld_pruned --genome --out ${SPECIES} --allow-extra-chr

## plot IBD
echo "pdf(\"${SPECIES}.IBD.pdf\",height=5,width=5)" > ${SPECIES}.ibd.plot
echo "d <- read.table(\"${SPECIES}.genome\",header=T)" >> ${SPECIES}.ibd.plot
echo "plot(d\$PI_HAT, d\$RATIO, cex = 0.5, xlab = \"PI_HAT\", ylab= \"RATIO\")" >> ${SPECIES}.ibd.plot 
echo "dev.off()" >> ${SPECIES}.ibd.plot
Rscript ${SPECIES}.ibd.plot

### Then, ADMIXTURE for k = 2 to 5 
### replace chromosome names in the bim file 
awk '{ print int($1), $2, $3, $4, $5, $6 }' < ${SPECIES}.ld_pruned.bim > ${SPECIES}.ld_pruned.bim.tmp
mv ${SPECIES}.ld_pruned.bim.tmp ${SPECIES}.ld_pruned.bim
for K in {2..5}
do
	admixture --cv ${SPECIES}.ld_pruned.bed $K > ${SPECIES}.${K}.admix.log 2> ${SPECIES}.${K}.admix.err
done

### collect admixture CVs for plot
cat ${SPECIES}.*.admix.log | grep CV | perl -pi -e 's/.+=//' | perl -pi -e 's/\): /\t/' > ${SPECIES}.CV

### now barplots for admixture output
echo "x <- read.table(\"${SPECIES}.CV\")" > ${SPECIES}.admixture.plot.r
echo "pdf(\"${SPECIES}.admix.pdf\",height=10,width=5)" >> ${SPECIES}.admixture.plot.r
echo "par(mfrow=c(5,1),mar=c(4,4,2,2))" >> ${SPECIES}.admixture.plot.r
echo "plot(x\$V1,x\$V2,xlab=\"K\",ylab=\"CV\")" >> ${SPECIES}.admixture.plot.r
echo "Q2 <- as.matrix(read.table(\"${SPECIES}.ld_pruned.2.Q\"))" >> ${SPECIES}.admixture.plot.r
echo "Q3 <- as.matrix(read.table(\"${SPECIES}.ld_pruned.3.Q\"))" >> ${SPECIES}.admixture.plot.r
echo "Q4 <- as.matrix(read.table(\"${SPECIES}.ld_pruned.4.Q\"))" >> ${SPECIES}.admixture.plot.r
echo "Q5 <- as.matrix(read.table(\"${SPECIES}.ld_pruned.5.Q\"))" >> ${SPECIES}.admixture.plot.r
echo "barplot(t(Q2),col=rainbow(2))" >> ${SPECIES}.admixture.plot.r
echo "barplot(t(Q3),col=rainbow(3))" >>	${SPECIES}.admixture.plot.r
echo "barplot(t(Q4),col=rainbow(4))" >>	${SPECIES}.admixture.plot.r
echo "barplot(t(Q5),col=rainbow(5))" >> ${SPECIES}.admixture.plot.r
echo "dev.off()" >> ${SPECIES}.admixture.plot.r

### run it
Rscript ${SPECIES}.admixture.plot.r
