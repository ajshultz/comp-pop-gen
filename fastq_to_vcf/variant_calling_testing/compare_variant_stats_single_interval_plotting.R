require(tidyverse)

setwd("~/Git_Repositories/comp-pop-gen/fastq_to_vcf/variant_calling_testing/")

tgut <- read_delim("_EXPT_COMPARE_STATS/Tguttata_int_1_stats_GATK_filters.txt",delim="\t")

tgut %>%
  mutate(pipeline = case_when(pipeline == "LC" ~ "low coverage",
                              pipeline == "HC" ~ "high coverage")) %>%
  mutate(coverage=paste0(coverage,"x")) %>%
  unite(coverage,pipeline,sep=" ",remove=FALSE,col="pipeline_coverage") %>%
  ggplot(aes(factor(pipeline_coverage,levels=c("2x low coverage","4x low coverage","15x low coverage","15x high coverage","30x low coverage","30x high coverage")),count)) +
  geom_bar(stat="identity",fill="#332288") +
  facet_wrap(~TYPE,nrow=1) +
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1)) +
  xlab("dataset")
ggsave("Tguttata_variant_counts.pdf",width=6,height=3)

tgut %>%
  mutate(pipeline = case_when(pipeline == "LC" ~ "low coverage",
                              pipeline == "HC" ~ "high coverage")) %>%
  mutate(coverage=paste0(coverage,"x")) %>%
  unite(coverage,pipeline,sep=" ",remove=FALSE,col="pipeline_coverage") %>%
  ggplot(aes(factor(pipeline_coverage,levels=c("2x low coverage","4x low coverage","15x low coverage","15x high coverage","30x low coverage","30x high coverage")),mean_het)) +
  geom_bar(stat="identity",fill="#332288") +
  facet_wrap(~TYPE,nrow=1) +
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1)) +
  xlab("dataset") +
  ylab("heterozygosity")
ggsave("Tguttata_heterozygosity.pdf",width=6,height=3)
