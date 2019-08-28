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

## cyvcf2.py output manipulations will take place here before 06_ 
