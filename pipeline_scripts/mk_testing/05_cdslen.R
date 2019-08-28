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

## as of 16:30 Aug 28, below is untested

corCor <- read.delim("corCor.effects.txt", sep = "\t", header = F, stringsAsFactors = F, quote = "") %>%
  as_tibble()
colnames(corCor) <- c("chr", "end.pos", "effect")
corCor <- corCor %>%
  mutate(start.pos = end.pos - 1) %>%
  select(chr, start.pos, end.pos, effect) 
write.table(corCor, "corCor.vcfann.bed", sep = "\t", row.names = F, quote = F)
