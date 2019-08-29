## Use of missingness data to determine if there are individuals with high relative missingness that need to be removed

install.packages("tidyverse")
libary(tidyverse)
corCor <- read.delim("/scratch/swuitchik/CompPopGen/missingness_ingroups/_Ccornix_all_all_missingness_info.txt")
corCor.clean <- select(corCor, -c(INTERVAL, N_GENOTYPES_FILTERED, F_MISS)) %>% 
  group_by(INDV) %>%
  summarise_each(funs(sum)) %>%
  mutate(missing = N_MISS/N_DATA)
corCor.threshold <- 3*median(corCor.clean$missing)
corCor.remove <- corCor.clean %>% 
  filter(missing >= corCor.threshold)
write.csv(corCor.remove %>% select(INDV), "corCor.remove.indv")

## if there are individuals to remove, use: vcftools --gzvcf corCor.clean.vcf.gz --remove-indv corCor.remove.indv --recode --recode-INFO-all --out corCor.clean2 
## use the VCF with high missingness individuals removed as input for subsequent step
