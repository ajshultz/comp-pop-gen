install.packages("tidyverse")
libary(tidyverse)
corCor <- read.delim("~/Desktop/PDF/CompPopGen/_Ccornix_all_all_missingness_info.txt")
corCor.clean <- select(corCor, -c(INTERVAL, N_GENOTYPES_FILTERED, F_MISS)) %>% 
  group_by(INDV) %>%
  summarise_each(funs(sum)) %>%
  mutate(Missing = N_MISS/N_DATA)
corCor.threshold <- 3*median(corCor.clean$Missing)
corCor.remove <- corCor.clean %>% 
  filter(Missing >= corCor.threshold)
write.csv(corCor.remove %>% select(INDV), "corCor_remove.csv")
