library(tidyverse)
library(SNPRelate)

##SNPRelate
setwd("~/Projects/popgen/softsweep/pca_test/")
vcf.fn <- "Mgallopavo.clean.recode.vcf.gz"
snpgdsVCF2GDS(vcf.fn, "Mgallopavo.gds", method="biallelic.only")
snpgdsSummary("Mgallopavo.gds")
genofile <- snpgdsOpen("Mgallopavo.gds")
snpset <- snpgdsLDpruning(genofile, autosome.only = FALSE, ld.threshold = 0.2)
snpset.id <- unlist(snpset)
pca <- snpgdsPCA(genofile,num.thread = 8, snp.id=snpset.id, autosome.only = FALSE, maf = 0.05, missing.rate = 0.1)

#pca plot
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
plot(tab$EV2, tab$EV1, xlab="", ylab="", pch=16, col="salmon", bty="l", cex=1.5, cex.axis=1.5, las=1)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

#identity by state MDS plot
ibs <- snpgdsIBS(genofile, num.thread=2, autosome.only = FALSE, missing.rate = 0)
image(ibs$ibs, col=terrain.colors(16))
loc <- cmdscale(1 - ibs$ibs, k = 2)
x <- loc[, 1]; y <- loc[, 2]
plot(x, y, xlab = "", ylab = "",
     main = "Multidimensional Scaling Analysis (IBS)")

