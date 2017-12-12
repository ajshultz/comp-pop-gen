library(tidyverse)


setwd("/n/holylfs/LABS/informatics/ashultz/CompPopGen")

args <- commandArgs(TRUE)
species <- args[1]

test.fn <- "~/Dropbox/Informatics/CompPopGen/Agambiae.depth_hist"
species <- "Agambiae"

depthhist.fn <- paste("/n/holylfs/LABS/informatics/tsackton/popgen/softsweep/stats/",species,".depth_hist",sep = "")
depthhist <- read_delim(test.fn,delim="\t",col_names = c("depth","count"))

nsites <- depthhist %>%
  summarise(sites=sum(count)) %>% pull

fewsites <- nsites*0.00001

perc_0 <- depthhist %>%
  mutate(perc=count/nsites*100) %>%
  mutate(perc=round(perc,1)) %>%
  filter(depth==0) %>%
  pull(perc)

max_depth_plot <- depthhist %>%
  filter(count<fewsites) %>%
  head(n=1) %>%
  pull(depth)

ggplot(depthhist,aes(depth,count)) +
  geom_col() +
  xlim(c(1,max_depth_plot)) +
  labs(title=paste(species,", ",nsites," sites, ",perc_0,"% missing data",sep=""))
ggsave(filename=paste0("depth_plots/",species,"restricted_dist.pdf"),device="pdf")

ggplot(depthhist,aes(depth,count)) +
  geom_col() +
  labs(title=paste(species,", ",nsites," sites, ",perc_0,"% missing data",sep=""))
ggsave(filename=paste0("depth_plots/",species,"entire_dist.pdf"),device="pdf")


