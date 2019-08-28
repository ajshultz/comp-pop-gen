rm(list=ls(all=TRUE))
cds <- read.delim("corCor.onlyCDS.genes.bed", sep = "\t", header = F, stringsAsFactors = F, quote = "") %>%
  as_tibble()
colnames(cds) <- c("chr", "start", "end", "gene")
cds.temp <- cds %>%
  mutate(cds.temp = end - start)
cds.len <- cds.temp %>%
  group_by(gene) %>% 
  summarise(cds.len = sum(cds.temp))
write.table(cds.len, "cdslen.txt", sep = "\t", row.names = F, quote = F)

corCor <- read.delim("corCor.ann.bed", header = F, sep = "\t", stringsAsFactors = F, quote = "") %>% 
  as_tibble()
colnames(corCor) <- c("#chr","start.pos","end.pos","id","qual","ref","alt","filter","info","i","ii","iii","iv","v","vi")
corCor <- corCor %>% select(-c(id, qual, ref, alt, filter, i, ii, iii, iv, v, vi))
corCor.clean <- data %>% separate(info, into = c(NA, "effect"), sep = "([\\|])")
write.table(corCor.clean, "corCor.vcfann.bed", sep = "\t", row.names = F, quote = F)

## NB: the cleaning of the annotated vcf will be streamlined with cyvcf2 by Aug 29
