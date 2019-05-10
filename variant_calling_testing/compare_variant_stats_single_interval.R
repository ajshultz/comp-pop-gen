require(tidyverse)

setwd("/Users/tim/Documents/Science/presentations/meetings/Cod Kickoff 2018")

species <- "Mmulatta"
species <- "Tguttata"
species <- "Herato"
species <- "Pacticauda"

LC_2X <- read_delim(paste0(species,"_2X_LC.txt.gz"),delim="\t") %>%
  rename(HOM_REF=`HOM-REF`,HOM_VAR=`HOM-VAR`,NO_CALL=`NO-CALL`) %>%
  mutate(EXP="LC_2X")
LC_4X <- read_delim(paste0(species,"_4X_LC.txt.gz"),delim="\t") %>%
  rename(HOM_REF=`HOM-REF`,HOM_VAR=`HOM-VAR`,NO_CALL=`NO-CALL`) %>%
  mutate(EXP="LC_4X")
LC_15X <- read_delim(paste0(species,"_15X_LC.txt.gz"),delim="\t") %>%
  rename(HOM_REF=`HOM-REF`,HOM_VAR=`HOM-VAR`,NO_CALL=`NO-CALL`) %>%
  mutate(EXP="LC_15X")
HC_15X <- read_delim(paste0(species,"_15X_HC.txt.gz"),delim="\t") %>%
  rename(HOM_REF=`HOM-REF`,HOM_VAR=`HOM-VAR`,NO_CALL=`NO-CALL`) %>%
  mutate(EXP="HC_15X")
LC_30X <- read_delim(paste0(species,"_30X_LC.txt.gz"),delim="\t") %>%
  rename(HOM_REF=`HOM-REF`,HOM_VAR=`HOM-VAR`,NO_CALL=`NO-CALL`) %>%
  mutate(EXP="LC_30X")
HC_30X <- read_delim(paste0(species,"_30X_HC.txt.gz"),delim="\t") %>%
  rename(HOM_REF=`HOM-REF`,HOM_VAR=`HOM-VAR`,NO_CALL=`NO-CALL`) %>%
  mutate(EXP="HC_30X")

all_res <- bind_rows(LC_2X,LC_4X,LC_15X,HC_15X,LC_30X,HC_30X)

#Get means of statistics with no filtering
all_res %>%
  group_by(EXP,TYPE) %>%
  summarize(count=n(),
            mean_het=mean(HET),
            mean_hom_ref=mean(HOM_REF),
            mean_hom_var=mean(HOM_VAR),
            mean_no_call=mean(NO_CALL),
            mean_QD=mean(QD,na.rm = T),
            mean_MQ=mean(MQ,na.rm = T),
            mean_FS=mean(FS,na.rm = T),
            mean_SOR=mean(SOR,na.rm = T)) %>%
  separate(EXP,into=c("pipeline","coverage")) %>%
  separate(coverage,into=c("coverage"),sep = "X",extra="drop") %>%
  mutate(coverage=as.numeric(coverage)) %>%
  arrange(desc(TYPE),coverage) %>%
  write_delim(paste0("_EXPT_COMPARE_STATS/",species,"_int_",interval,"_stats_no_filter.txt"),delim="\t")

#Get means of statistics with standard filtering
all_res_filtered_SNPs <- all_res %>%
  group_by(EXP,TYPE) %>%
  filter(TYPE=="SNP", !(QD < 2 | FS > 60.0 | MQ < 40 | SOR > 3.0 | MQRankSum < -12.5 | ReadPosRankSum < -8)) %>%
  summarize(count=n(),
            mean_het=mean(HET),
            mean_hom_ref=mean(HOM_REF),
            mean_hom_var=mean(HOM_VAR),
            mean_no_call=mean(NO_CALL),
            mean_QD=mean(QD,na.rm = T),
            mean_MQ=mean(MQ,na.rm = T),
            mean_FS=mean(FS,na.rm = T),
            mean_SOR=mean(SOR,na.rm = T)) %>%
  separate(EXP,into=c("pipeline","coverage")) %>%
  separate(coverage,into=c("coverage"),sep = "X",extra="drop") %>%
  mutate(coverage=as.numeric(coverage)) %>%
  arrange(desc(TYPE),coverage)

all_res_filtered_INDELs <- all_res %>%
  group_by(EXP,TYPE) %>%
  filter(TYPE=="INDEL", !(QD < 2 | FS > 200.0 | SOR > 10.0 | ReadPosRankSum > 20)) %>%
  summarize(count=n(),
            mean_het=mean(HET),
            mean_hom_ref=mean(HOM_REF),
            mean_hom_var=mean(HOM_VAR),
            mean_no_call=mean(NO_CALL),
            mean_QD=mean(QD,na.rm = T),
            mean_MQ=mean(MQ,na.rm = T),
            mean_FS=mean(FS,na.rm = T),
            mean_SOR=mean(SOR,na.rm = T)) %>%
  separate(EXP,into=c("pipeline","coverage")) %>%
  separate(coverage,into=c("coverage"),sep = "X",extra="drop") %>%
  mutate(coverage=as.numeric(coverage)) %>%
  arrange(desc(TYPE),coverage)

all_res_filtered_both <- bind_rows(all_res_filtered_SNPs,all_res_filtered_INDELs)

all_res_filtered_both %>%
  write_delim(paste0("_EXPT_COMPARE_STATS/",species,"_int_",interval,"_stats_GATK_filters.txt"),delim="\t")

all_res %>%
  filter(TYPE=="SNP") %>%
  separate(EXP,into=c("pipeline","coverage")) %>%
  separate(coverage,into=c("coverage"),sep = "X",extra="drop") %>%
  mutate(coverage=as.numeric(coverage)) %>%
  ggplot(aes(QD)) +
  geom_density() +
  facet_grid(coverage~pipeline)
ggsave(paste0("_EXPT_COMPARE_STATS/",species,"_int_",interval,"_SNP_QD_plots.pdf"))

all_res %>%
  filter(TYPE=="SNP") %>%
  separate(EXP,into=c("pipeline","coverage")) %>%
  separate(coverage,into=c("coverage"),sep = "X",extra="drop") %>%
  mutate(coverage=as.numeric(coverage)) %>%
  ggplot(aes(MQ)) +
  geom_density() +
  xlim(c(20,80)) +
  facet_grid(coverage~pipeline)
ggsave(paste0("_EXPT_COMPARE_STATS/",species,"_int_",interval,"_SNP_MQ_plots.pdf"))

all_res %>%
  filter(TYPE=="SNP") %>%
  separate(EXP,into=c("pipeline","coverage")) %>%
  separate(coverage,into=c("coverage"),sep = "X",extra="drop") %>%
  mutate(coverage=as.numeric(coverage)) %>%
  ggplot(aes(FS)) +
  geom_density() +
  facet_grid(coverage~pipeline)
ggsave(paste0("_EXPT_COMPARE_STATS/",species,"_int_",interval,"_SNP_FS_plots.pdf"))

all_res %>%
  filter(TYPE=="SNP") %>%
  separate(EXP,into=c("pipeline","coverage")) %>%
  separate(coverage,into=c("coverage"),sep = "X",extra="drop") %>%
  mutate(coverage=as.numeric(coverage)) %>%
  ggplot(aes(SOR)) +
  geom_density() +
  facet_grid(coverage~pipeline)
ggsave(paste0("_EXPT_COMPARE_STATS/",species,"_int_",interval,"_SNP_SOR_plots.pdf"))

all_res %>%
  filter(TYPE=="SNP") %>%
  separate(EXP,into=c("pipeline","coverage")) %>%
  separate(coverage,into=c("coverage"),sep = "X",extra="drop") %>%
  mutate(coverage=as.numeric(coverage)) %>%
  ggplot(aes(MQRankSum)) +
  geom_density() +
  facet_grid(coverage~pipeline)
ggsave(paste0("_EXPT_COMPARE_STATS/",species,"_int_",interval,"_SNP_MQRankSum_plots.pdf"))

all_res %>%
  filter(TYPE=="SNP") %>%
  separate(EXP,into=c("pipeline","coverage")) %>%
  separate(coverage,into=c("coverage"),sep = "X",extra="drop") %>%
  mutate(coverage=as.numeric(coverage)) %>%
  ggplot(aes(ReadPosRankSum)) +
  geom_density() +
  facet_grid(coverage~pipeline)
ggsave(paste0("_EXPT_COMPARE_STATS/",species,"_int_",interval,"_SNP_ReadPosRankSum_plots.pdf"))



