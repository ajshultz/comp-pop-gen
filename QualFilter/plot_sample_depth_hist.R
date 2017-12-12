library(tidyverse)


setwd("/n/holylfs/LABS/informatics/ashultz/CompPopGen")

args <- commandArgs(TRUE)
species <- args[1]

#test.fn <- "~/Dropbox/Informatics/CompPopGen/Agambiae.samp_depth"
#species <- "Agambiae"
#depthhist <- read_delim(test.fn,delim="\t",col_names = c("depth","count"))

depthhist.fn <- paste("/n/holylfs/LABS/informatics/tsackton/popgen/softsweep/stats/",species,".depth_hist",sep = "")
depthhist <- read_delim(depthhist.fn,delim="\t",col_names = c("depth","count"))

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
ggsave(filename=paste0("depth_plots/",species,"_mean_sample_depth.pdf"),device="pdf")



