library(tidyverse)
library(SNPRelate)

setwd("/n/holylfs/LABS/informatics/ashultz/CompPopGen")

args <- commandArgs(TRUE)
species <- args[1]

##SNPRelate
vcf.fn <- paste("/n/holylfs/LABS/informatics/tsackton/popgen/softsweep/clean_vcfs/",species,".clean.recode.vcf.gz",sep = "")
gds.fn <- paste("gds/",species,".gds",sep = "")
snpgdsVCF2GDS(vcf.fn, gds.fn, method="biallelic.only")
snpgdsSummary(gds.fn)
genofile <- snpgdsOpen(gds.fn)
snpset <- snpgdsLDpruning(genofile, autosome.only = FALSE, ld.threshold = 0.2)
snpset.id <- unlist(snpset)
pca <- snpgdsPCA(genofile,num.thread = 8, snp.id=snpset.id, autosome.only = FALSE, maf = 0.05, missing.rate = 0.1)

#pca plot
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
pdf(file=paste("pca_plots/",species,"_pc1pc2.pdf",sep=""))
pc.percent <- pca$varprop*100 %>%
  round(.,2)
plot(tab$EV2, tab$EV1, xlab=paste("PC2 ",pc.percent[2],"%",sep=""), ylab=paste("PC1 ",pc.percent[1],"%",sep=""), pch=16, col="salmon", bty="l", cex=1.5, cex.axis=1.5, las=1)
dev.off()

#identity by state MDS plot
ibs <- snpgdsIBS(genofile, num.thread=2, autosome.only = FALSE, missing.rate = 0)
image(ibs$ibs, col=terrain.colors(16))
loc <- cmdscale(1 - ibs$ibs, k = 2)
x <- loc[, 1]; y <- loc[, 2]
pdf(file=paste("pca_plots/",species,"_ibs.pdf",sep=""))
plot(x, y, xlab = "", ylab = "",
     main = "Multidimensional Scaling Analysis (IBS)")
dev.off()

