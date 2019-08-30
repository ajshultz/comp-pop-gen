## Final table creation

cds <- read.delim("corCor.onlyCDS.genes.bed", sep = "\t", header = F, stringsAsFactors = F, quote = "") %>%
  as_tibble()
colnames(cds) <- c("chr", "start", "end", "gene")
cds.temp <- cds %>%
  mutate(cds.temp = end - start)
cds.len <- cds.temp %>%
  group_by(gene) %>% 
  summarise(cds.len = sum(cds.temp))
write.table(cds.len, "cdslen.txt", sep = "\t", row.names = F, quote = F)

call <- read.delim("corCor.callable.bed", sep = "\t", header = F, stringsAsFactors = F, quote = "") %>% as_tibble()
colnames(call) <- c("chr", "start", "end", "gene")
call.temp <- call %>% mutate(call.temp = end - start)
call.length <- call.temp %>% group_by(gene) %>% summarise(call.len = sum(call.temp))
write.table(call.length, "call.txt", sep = "\t", row.names = F, quote = F)

data <- read.delim("corCor.final.clean.bed", sep = "\t", header = F, stringsAsFactors = F, quote = "") %>% as_tibble()
colnames(data) <- c("chr", "start", "end", "gene", "effect")
mis <- data %>% 
  group_by(gene) %>% 
  tally(effect ==  "missense_variant")
colnames(mis) <- c("gene", "mis")
syn <- data %>%
  group_by(gene) %>%
  tally(effect == "synonymous_variant")
colnames(syn) <- c("gene", "syn")

table <- left_join(mis, syn, by = "gene")
tab2 <- left_join(table, cds.len, by = "gene")
final <- left_join(tab2, call, by = "gene") %>% 
  select(gene, cds.len, call.length, syn, mis)

# qc: make sure no ratios > 1
check <- final %>% 
  mutate(check = call.len/cds.len)
check[!complete.cases(check),]
max(check$check)
write.table(final, "corCor.table.txt", sep = "\t", row.names = F, quote = F)
# write.table(check, "corCor.table_check.txt", sep = "\t", row.names = F, quote = F)
