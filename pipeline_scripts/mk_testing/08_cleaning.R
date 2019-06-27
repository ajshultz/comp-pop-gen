library(tidyverse)
consensus <- read.table("corCor.consensus.bed", header = F, sep = "\t", stringsAsFactors = F, quote = "") %>%
  as_tibble()
colnames(consensus, do.NULL = T, prefix = "col")
colnames(consensus) <- c("chr", "start.pos", "end.pos","info")
