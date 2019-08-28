call <- read.delim("corCor.callable.bed", sep = "\t", header = F, stringsAsFactors = F, quote = "") %>% as_tibble()
colnames(call) <- c("chr", "start", "end", "gene")
call.temp <- call %>% mutate(call.temp = end - start)
call.len <- call.temp %>% group_by(gene) %>% summarise(call.len = sum(call.temp))
write.table(call.len, "call.txt", sep = "\t", row.names = F, quote = F)

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
cds <- read.delim("cdslen.txt", sep = "\t", stringsAsFactors = F, quote = "") %>% as_tibble()
tab2 <- left_join(table, cds, by = "gene")
call <- read.delim("call.txt", sep = "\t", stringsAsFactors = F, quote = "") %>% as_tibble()
final <- left_join(tab2, call, by = "gene") %>% select(gene, cds.len, call.len, syn, mis)
check <- final %>% mutate(check = call.len/cds.len)
check[!complete.cases(check),]
max(check$check)
write.table(final, "corCor.table.txt", sep = "\t", row.names = F, quote = F)
write.table(check, "corCor.table_check.txt", sep = "\t", row.names = F, quote = F)
