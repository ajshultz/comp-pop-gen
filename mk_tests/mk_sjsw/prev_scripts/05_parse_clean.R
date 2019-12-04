## Clean effects file and convert to a BED format

corCor <- read.delim("corCor.effects.txt", sep = "\t", header = F, stringsAsFactors = F, quote = "") %>%
  as_tibble()
colnames(corCor) <- c("#chr", "end.pos", "effect")
corCor <- corCor %>%
  mutate(start.pos = end.pos - 1) %>%
  select('#chr', start.pos, end.pos, effect)
write.table(corCor, "corCor.vcfann.bed", sep = "\t", row.names = F, quote = F)
