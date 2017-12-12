library(tidyverse)


setwd("/n/holylfs/LABS/informatics/ashultz/CompPopGen")

args <- commandArgs(TRUE)
species <- args[1]

#test.fn <- "~/Dropbox/Informatics/CompPopGen/Agambiae.samp_depth"
#species <- "Agambiae"
#depthhist <- read_delim(test.fn,delim="\t",col_names = c("depth","count"))

depthhist.fn <- paste("/n/holylfs/LABS/informatics/tsackton/popgen/softsweep/stats/",species,".samp_depth",sep = "")
depthhist <- read_delim(depthhist.fn,delim="\t",col_names = c("depth","count"))

nsites <- depthhist %>%
  summarise(sites=sum(count)) %>% pull

perc_0 <- depthhist %>%
  mutate(perc=count/nsites*100) %>%
  mutate(perc=round(perc,1)) %>%
  filter(depth==0) %>%
  pull(perc)


depthhist %>%
  mutate(ints=round(depth,0)) %>%
  group_by(ints) %>%
  summarize(max_depth = mean(count)) %>%
  ggplot(aes(ints,max_depth)) +
  geom_col() +
  xlim(c(1,100)) +
  xlab("mean depth") +
  ylab("counts") +
  labs(title=paste(species,", ",nsites," sites, ",perc_0,"% missing data",sep=""))
ggsave(filename=paste0("sample_depth_plots/",species,"_mean_sample_depth_100.pdf"),device="pdf")


depthhist %>%
  mutate(ints=round(depth,0)) %>%
  group_by(ints) %>%
  summarize(max_depth = mean(count)) %>%
  ggplot(aes(ints,max_depth)) +
  geom_col() +
  xlim(c(1,25)) +
  xlab("mean depth") +
  ylab("counts") +
  labs(title=paste(species,", ",nsites," sites, ",perc_0,"% missing data",sep=""))
ggsave(filename=paste0("sample_depth_plots/",species,"_mean_sample_depth_25.pdf"),device="pdf")

