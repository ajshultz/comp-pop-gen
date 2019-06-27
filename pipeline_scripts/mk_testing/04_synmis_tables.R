library(tidyverse)
library(tictoc)
library(beepr)
corCor <- read.delim("~/Desktop/PDF/CompPopGen/annotated/corCor.oneper.txt") %>% 
  as_tibble()
corCor.clean <- corCor %>% 
  rename(chrom = "CHROM", pos = "POS", ann.effect = "ANN....EFFECT", ann.geneid = "ANN....GENEID")
tic("mis")
mis <- corCor.clean %>% 
  group_by(ann.geneid) %>% 
  tally(ann.effect ==  "missense_variant")
toc(); beep("coin") 
mis.count <- sum(mis$n)
mis.table <- corCor.clean %>% 
  filter(ann.effect =="missense_variant") %>%
  group_by(ann.geneid) 
syn <- corCor.clean %>% 
  group_by(ann.geneid) %>% 
  tally(ann.effect ==  "synonymous_variant"); 
syn.count <- sum(syn$n)
syn.table <- corCor.clean %>% 
  filter(ann.effect =="synonymous_variant")
write.csv(syn.table, "corCor.syn.csv")
write.csv(mis.table, "corCor.mis.csv")
