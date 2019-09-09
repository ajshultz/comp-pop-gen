setwd("~/Dropbox/BirdImmuneGeneEvolution/PopGen/")

library(ggplot2)
library(gridExtra)
library(tidyverse)


#Calculate standard MK test, returns a list with the pvalue and alpha value
mk_test <- function(dn=dn,ds=ds,pn=pn,ps=ps){
  dnds_mat <- matrix(data=c(ds,dn,ps,pn),nrow=2,byrow = F)
  
  #Fisher's exact test to caluclate significance
  pval = fisher.test(dnds_mat)$p.value
  #Calculate alpha, or proportion of sites estimated to be under positive selection
  alpha = 1-((ds*pn)/(dn*ps))
  
  return(list(pval,alpha))
}

#Input a result object from snipre (as long as it contains snipre input columns), and return that object including traditional MK alpha values, p-values and FDR-corrected p-values.
mk_tibble_calc <- function(snipre_res_obj){
  snipre_res_obj_new <- snipre_res_obj %>%
    rowwise %>%
    mutate(mk_pval = mk_test(dn=FR,ds=FS,pn=PR,ps=PS)[[1]],
           alpha = mk_test(dn=FR,ds=FS,pn=PR,ps=PS)[[2]]) %>%
    ungroup %>%
    dplyr::mutate(mk_pval_fdr = p.adjust(mk_pval,method="BH"))
  return(snipre_res_obj_new)
}

#Given ML and Bayesian file locations and an abbreviation, will read in ML results and Bayesian results (assuming created with standard script in standard locations), select necessary columns of interest, conduct standard MK tests, and combine them based on gene symbol, producing a single tibble.
proc_ml_bayes_res <- function(sp_abbr,ml_file,bayes_file){
  
  sp_ml_res <- read_csv(ml_file)
  sp_bay_res <- read_csv(bayes_file)
  
  sp_bay_res_sim <- sp_bay_res %>%
    dplyr::select(geneID,BSnIPRE.class,BSnIPRE.est,BSnIPRE.lbound,BSnIPRE.ubound)
  
  sp_ml_res <- mk_tibble_calc(sp_ml_res)
  
  sp_combo_res <- sp_ml_res %>%
    mutate(species = sp_abbr) %>%
    dplyr::select(species,geneID,PR,FR,PS,FS,Tsil,Trepl,nout,npop,SnIPRE.class,SnIPRE.est,SnIPRE.lbound,SnIPRE.ubound,alpha,mk_pval,mk_pval_fdr) %>%
    left_join(sp_bay_res_sim) %>%
    dplyr::rename(gene_name = geneID)
  
  return(sp_combo_res)
}

#abbreviations of all species tested
sp_abbr <- list("Gvarius","Cjaponica","Tguttata","Pmajor","Egarzetta","Nnippon","Ccornix")

#list to hold all resulting mk tibbles
mk_res_list <- list()

#Run through all species and process ml and Bayesian result files
for (i in 1:length(sp_abbr)){
  mk_res_list[[i]] <- proc_ml_bayes_res(sp_abbr = sp_abbr[[i]],ml_file = paste0("TestSpecies_2019/",sp_abbr[[i]],"_results.csv"),bayes_file = paste0("TestSpecies_2019/SnIPRE_code_JAGS/",sp_abbr[[i]],"_bayesian_results.csv"))
  
}

#Cat all results into a single tibble
mk_res <- bind_rows(mk_res_list)


#Falbicollis with new machinery would only run for Bayesian file. Run that and bring into dataset.
fa_bay_res <- read_csv("TestSpecies_2019/SnIPRE_code_JAGS/Falbicollis_bayesian_results.csv")

fa_bay_res_proc <- fa_bay_res %>%
  mk_tibble_calc %>%
  mutate(species="Falbicollis",SnIPRE.class=NA,SnIPRE.est=NA,SnIPRE.lbound=NA,SnIPRE.ubound=NA) %>%
  dplyr::rename(gene_name = geneID) %>%
  dplyr::select(species,gene_name,PR,FR,PS,FS,Tsil,Trepl,nout,npop,SnIPRE.class,SnIPRE.est,SnIPRE.lbound,SnIPRE.ubound,alpha,mk_pval,mk_pval_fdr,BSnIPRE.class,BSnIPRE.est,BSnIPRE.lbound,BSnIPRE.ubound)

mk_res <- bind_rows(mk_res,fa_bay_res_proc)

#Load data for Falbicollis as calculated with previous VCF file last year
fa_ml_res <- read_csv("TestSpecies/fa_eb_results.csv")
fa_bay_res <- read_csv("TestSpecies/Falbicolis.Tguttata.maf.snipre.bayesianresults.csv")
fa_rna_to_gene <- read_delim("TestSpecies/all_info.summary_ficAlb",delim="\t",col_names = c("species",	"gene_id","gene_acc","gene_name","biotype","trans_id","trans_acc","trans_len","cds_len","prot_id","is_longest_prot"))

fa_bay_res_sim <- fa_bay_res %>%
  dplyr::select(geneID,BSnIPRE.class,BSnIPRE.est,BSnIPRE.lbound,BSnIPRE.ubound)

fa_ml_res <- mk_tibble_calc(fa_ml_res)

fa_combo_res <- fa_ml_res %>%
  mutate(species="Falbicollis_Previous") %>%
  dplyr::select(species,geneID,PR,FR,PS,FS,Tsil,Trepl,nout,npop,SnIPRE.class,SnIPRE.est,SnIPRE.lbound,SnIPRE.ubound,alpha,mk_pval,mk_pval_fdr) %>%
  left_join(fa_bay_res_sim) %>%
  dplyr::rename(trans_id = geneID)

#Create hog ID maps
hog_id_map <- read_delim("../new_hog_list.txt",delim="\t",col_names = c("hog","ensembl","gene_acc","species"))
hog_id_gg <- hog_id_map %>% filter(species=="galGal") %>% dplyr::select(hog,entrezgene=gene_acc)

#Get geneIDs, only keep longest transcripts (those used in alignment), and drop duplicate proteins, join to get chicken gene IDs and therefore hogs.
fa_combo_res_renamed <- fa_combo_res %>%
  left_join(fa_rna_to_gene,by = "trans_id") %>% filter(is_longest_prot=="Y") %>%
  distinct(gene_acc,.keep_all=TRUE) %>%
  left_join(hog_id_map,by="gene_acc") %>% left_join(hog_id_gg,by="hog") %>%
  distinct(hog,.keep_all=TRUE) %>%
  dplyr::select(species=species.x,gene_name,PR,FR,PS,FS,Tsil,Trepl,nout,npop,SnIPRE.class,SnIPRE.est,SnIPRE.lbound,SnIPRE.ubound,alpha,mk_pval,mk_pval_fdr,BSnIPRE.class,BSnIPRE.est,BSnIPRE.lbound,BSnIPRE.ubound)


mk_res <- bind_rows(mk_res,fa_combo_res_renamed)


mk_res_summary_table <- mk_res %>%
  group_by(species) %>%
  summarize(n_sp = mean(npop),
            n_out=mean(nout),
            n_loci_tested = n(),
            n_sig_mk = sum(mk_pval<0.05),
            n_sig_mk_fdr = sum(mk_pval_fdr<0.05),
            n_pos_snipre_ml = sum(SnIPRE.class == "pos"),
            n_pos_bsnipre_bayes = sum(BSnIPRE.class == "pos"),
            n_pos_both_snipre = sum(SnIPRE.class == "pos" & BSnIPRE.class == "pos"),
            mean_SnIPRE.est = mean(SnIPRE.est),
            mean_BSnIPRE.est = mean(BSnIPRE.est)) %>%
  mutate(prop_pos_both = n_pos_both_snipre/n_loci_tested)

write_csv(mk_res_summary_table,"TestSpecies_2019/mk_res_summary_table_all_sp.csv")

#Read in original polymorphism stats for visualization purposes

species <- "Tguttata"
species <- "Falbicollis"
sp_abbr_fa <- c(sp_abbr,"Falbicollis")

poly_data_list <- list()

for (i in 1:length(sp_abbr_fa)){
  poly_data_list[[i]] <- read_delim(paste0("TestSpecies_2019/",sp_abbr_fa[[i]],"_PolymorphismDivergenceStats_combined.txt"),delim="\t") %>%
    dplyr::select(-X10) %>%
    mutate(NonsynSynSites = CallableNsynSites/CallableSynSites) %>%
    mutate(NonsynSynSNPs = NonsynPolymorphic/SynPolymorphic) %>%
    mutate(species = sp_abbr_fa[[i]])
}

poly_data <- bind_rows(poly_data_list)

#Count the number of sites with at least one callable site, and at least one SNP
poly_data %>%
  mutate(any_callable_sites = if_else(CallableSites>0,TRUE,FALSE)) %>%
  mutate(any_snps_present = if_else(SynPolymorphic>0|NonsynPolymorphic>0|SynFixed>0|NonsynFixed>0, TRUE,FALSE)) %>%
  group_by(species) %>%
  count(any_callable_sites,any_snps_present) %>%
  write_csv("TestSpecies_2019/callable_sites_snps_stats.csv")

poly_data %>%
  ggplot(aes(log10(CallableSites),fill=species)) +
  geom_histogram(bins=100) +
  facet_wrap(~species) +
  theme_bw()
ggsave("TestSpecies_2019/CallableSites_Histogram.pdf",width=10,height = 8)

poly_data %>%
  filter(SynPolymorphic>0|NonsynPolymorphic>0|SynFixed>0|NonsynFixed>0,CallableSites>0) %>%
  ggplot(aes(NonsynSynSites,fill=species)) +
  geom_histogram(bins=100) +
  facet_wrap(~species) +
  theme_bw() +
  labs(title = "Ratio of Nonsyn to Syn Callable Sites")
ggsave("TestSpecies_2019/NonsynSynCallableSites_Histogram.pdf",width=10,height = 8)

poly_data %>%
  filter(SynPolymorphic>0|NonsynPolymorphic>0|SynFixed>0|NonsynFixed>0,CallableSites>0) %>%
  filter(NonsynSynSites < 5) %>%
  ggplot(aes(NonsynSynSites,fill=species)) +
  geom_histogram(bins=100) +
  facet_wrap(~species) +
  theme_bw() +
  labs(title = "Ratio of Nonsyn to Syn Callable Sites Filtered")
ggsave("TestSpecies_2019/NonsynSynCallableSites_LessThan5_Histogram.pdf",width=10,height = 8)
  

poly_data %>%
  filter(SynPolymorphic>0|NonsynPolymorphic>0|SynFixed>0|NonsynFixed>0,CallableSites>0) %>%
  ggplot(aes(NonsynSynSNPs,fill=species)) +
  geom_histogram(bins=100) +
  facet_wrap(~species) +
  theme_bw() +
  labs(title = "Ratio of Nonsyn to Syn Poly SNPs")
ggsave("TestSpecies_2019/NonsynSynPolySNPs_Histogram.pdf",width=10,height = 8)

poly_data %>%
  filter(SynPolymorphic>0|NonsynPolymorphic>0|SynFixed>0|NonsynFixed>0,CallableSites>0) %>%
  filter(NonsynSynSNPs < 2) %>%
  ggplot(aes(NonsynSynSNPs,fill=species)) +
  geom_histogram(bins=100) +
  facet_wrap(~species) +
  theme_bw() +
  labs(title = "Ratio of Nonsyn to Syn Poly SNPs")
ggsave("TestSpecies_2019/NonsynSynPolySNPs_LessThan2_Histogram.pdf",width=10,height = 8)


####################
#Load comparative results
comp_res <- load("../02_output_annotated_data/all_res_zf_hs.Rdat")

#ZF gene ID table:
zf_gene_ids <- read_delim("~/Dropbox/HFWGReseq/GFF_Functional_Annotation/ZF_GeneID_Table.txt",delim="\t")
zf_gene_ids <- zf_gene_ids %>%
  dplyr::select(GeneID,Symbol,Aliases,description,chromosome)

all_res_gene_zf_hs <- all_res_gene_zf_hs %>%
  mutate(all_sel=if_else(FDRPval_busted<0.05 & FDRPval_m1m2<0.05 & FDRPval_m2m2a<0.05 & FDRPval_m7m8<0.05 & FDRPval_m8m8a<0.05,1,0),
         gene_name_both = coalesce(external_gene_name.x,external_gene_name.y))


all_sel_names <- all_res_gene_zf_hs %>%
  dplyr::filter(all_sel == 1) %>%
  dplyr::select(gene_name_both) %>%
  distinct(gene_name_both) %>%
  pull(gene_name_both)
all_names <- all_res_gene_zf_hs %>%
  dplyr::select(gene_name_both) %>%
  distinct(gene_name_both) %>%
  pull(gene_name_both)

mk_res_comp <- mk_res %>%
  mutate(sel_birds = if_else(gene_name %in% all_names,true = "tested",false="not_tested")) %>%
  mutate(sel_birds = if_else(gene_name %in% all_sel_names,"sig",sel_birds))
  
mk_res_comp %>%
  group_by(species) %>%
  count(sel_birds,BSnIPRE.class) %>%
  spread(BSnIPRE.class,n) %>%
  mutate(prop_sig_sel = pos/(sum(neut,pos,neg,na.rm = T))) %>%
  filter(sel_birds != "not_tested") %>%
  select(species,sel_birds,prop_sig_sel) %>%
  spread(sel_birds,prop_sig_sel) %>%
  mutate(enrichment = sig/tested)

mk_res_comp %>%
  mutate(mk_sig = if_else(SnIPRE.class == "pos" | BSnIPRE.class == "pos",1,0)) %>%
  group_by(species) %>%
  count(sel_birds,mk_sig) %>%
  filter(mk_sig == 1) %>%
  spread(sel_birds,n) %>%
  mutate(prop_sig_selbirds = (sig/sum(sig,tested,na.rm=T)))

mk_res_comp %>%
  mutate(mk_sig = if_else(SnIPRE.class == "pos" | BSnIPRE.class == "pos",1,0)) %>%
  group_by(gene_name) %>%
  summarize(n_sig_sp = sum(mk_sig)) %>%
  arrange(desc(n_sig_sp)) %>%
  mutate(sel_birds = if_else(gene_name %in% all_names,true = "tested",false="not_tested")) %>%
  mutate(sel_birds = if_else(gene_name %in% all_sel_names,"sig",sel_birds)) %>%
  filter(n_sig_sp >= 3)


#Pathway enrichment
library(DOSE)
library(clusterProfiler)
library(biomaRt)

mart <- useMart(biomart="ensembl",dataset="hsapiens_gene_ensembl")
hs_trans_table <- getBM(attributes =c("hgnc_symbol","entrezgene"), mart=mart)

#Pull all entrezgene IDs (human) for genes postiviely selected in 2 or more species
sig_2 <- mk_res_comp %>%
  mutate(mk_sig = if_else(SnIPRE.class == "pos" | BSnIPRE.class == "pos",1,0)) %>%
  group_by(gene_name) %>%
  summarize(n_sig_sp = sum(mk_sig)) %>%
  arrange(desc(n_sig_sp)) %>%
  mutate(sel_birds = if_else(gene_name %in% all_names,true = "tested",false="not_tested")) %>%
  mutate(sel_birds = if_else(gene_name %in% all_sel_names,"sig",sel_birds)) %>%
  filter(n_sig_sp >= 2) %>%
  left_join(hs_trans_table,by=c("gene_name" = "hgnc_symbol")) %>%
  pull(entrezgene)

#Pull all entrezgenes (human) tested in at least 3 species
all_tested <- mk_res_comp %>%
  group_by(gene_name) %>%
  summarize(n_sp = n()) %>%
  filter(n_sp >= 2) %>%
  left_join(hs_trans_table,by=c("gene_name" = "hgnc_symbol")) %>%
  pull(entrezgene)

#Set univeres of genes to those in results for all tests
pos2_genes_k <- enrichKEGG(sig_2,organism="hsa",pvalueCutoff=1,pAdjustMethod="BH",qvalueCutoff=.2,universe=as.character(all_tested),keyType="ncbi-geneid")

as.data.frame(pos2_genes_k)

cnetplot(pos2_genes_k,showCategory = 10)
dotplot(pos2_genes_k)
summary(pos2_genes_k)









#Combine tg and fa results by gene_name, combine with zf gene id table, combine with comp_res by gene name
fa_combo_res_renamed_simple <- fa_combo_res_renamed %>%
  dplyr::select(gene_name,fa_SnIPRE.class=SnIPRE.class,fa_BSnIPRE.class=BSnIPRE.class,fa_SnIPRE.est=SnIPRE.est,entrezgene,ensembl)

gv_combo_res_simple <- gv_combo_res %>%
  dplyr::select(gene_name,gv_SnIPRE.class=SnIPRE.class,gv_BSnIPRE.class=BSnIPRE.class,gv_SnIPRE.est=SnIPRE.est)

both_res <- tg_combo_res %>%
  dplyr::select(gene_name,tg_SnIPRE.class=SnIPRE.class,tg_BSnIPRE.class=BSnIPRE.class,tg_SnIPRE.est=SnIPRE.est) %>%
  full_join(fa_combo_res_renamed_simple, by = c("gene_name")) %>%
  full_join(gv_combo_res_simple,by=c("gene_name")) %>%
  left_join(zf_gene_ids,by=c("gene_name" = "Symbol")) %>%
  full_join(all_res_gene_zf_hs,by=c("gene_name" = "external_gene_name.x"))

both_res %>%
  filter(!is.na(tg_BSnIPRE.class),!is.na(fa_BSnIPRE.class)) %>%
  group_by(tg_BSnIPRE.class,fa_BSnIPRE.class,gv_BSnIPRE.class) %>%
  count()

both_res %>%
  filter(tg_BSnIPRE.class == "pos" | tg_SnIPRE.class == "pos", fa_BSnIPRE.class == "pos" | fa_BSnIPRE.class == "pos", gv_BSnIPRE.class == "pos" | gv_SnIPRE.class == "pos", all_sel ==1)




#Pathway enrichment
library(DOSE)
library(clusterProfiler)

fa_genes_signames <- fa_combo_res_renamed_simple %>%
  filter(fa_BSnIPRE.class == "pos" | fa_SnIPRE.class == "pos") %>%
  filter(!is.na(entrezgene)) %>%
  pull(entrezgene)
all_genes_fa <- fa_combo_res_renamed_simple %>%
  filter(!is.na(entrezgene),!is.na(fa_SnIPRE.class)) %>%
  pull(entrezgene)

#Set univeres of genes to those in results for all tests
fa_genes_k <- enrichKEGG(fa_genes_signames,organism="gga",pvalueCutoff=1,pAdjustMethod="BH",qvalueCutoff=1,universe=as.character(all_genes_fa),keyType="ncbi-geneid")

cnetplot(fa_genes_k)
dotplot(fa_genes_k)
summary(fa_genes_k)

as.data.frame(fa_genes_k)


fa_genes_signames <- both_res %>% filter(fa_BSnIPRE.class == "pos" | fa_SnIPRE.class == "pos") %>% pull(entrezgene.y)
fa_genes_signames <- fa_genes_signames[!is.na(fa_genes_signames)]
all_genes <- both_res %>% filter(!is.na(entrezgene.y),!is.na(fa_SnIPRE.class)) %>% pull(entrezgene.y)

#Set univeres of genes to those in results for all tests
fa_genes_k <- enrichKEGG(fa_genes_signames,organism="gga",pvalueCutoff=1,pAdjustMethod="BH",qvalueCutoff=1,universe=as.character(all_genes),keyType="ncbi-geneid")

cnetplot(fa_genes_k)
dotplot(fa_genes_k)
summary(fa_genes_k)
write_csv(as.data.frame(fa_genes_k),"TestSpecies/fa_pathwayEnrichment.csv")


#ZF from ZF data alone
tg_genes_signames <- tg_combo_res %>%
  left_join(zf_gene_ids,by=c("gene_name" = "Symbol")) %>%
  filter(BSnIPRE.class == "pos" | SnIPRE.class == "pos") %>%
  filter(!is.na(GeneID)) %>%
  pull(GeneID)
all_genes_tg <- tg_combo_res %>%
  left_join(zf_gene_ids,by=c("gene_name" = "Symbol")) %>%
  filter(!is.na(GeneID),!is.na(SnIPRE.class)) %>%
  pull(GeneID)

#Set univeres of genes to those in results for all tests
tg_genes_k <- enrichKEGG(tg_genes_signames,organism="tgu",pvalueCutoff=1,pAdjustMethod="BH",qvalueCutoff=1,universe=as.character(all_genes_tg),keyType="ncbi-geneid")

cnetplot(tg_genes_k)
dotplot(tg_genes_k)
summary(tg_genes_k)

write_csv(as.data.frame(tg_genes_k),"TestSpecies_2019/Tguttata_pathwayEnrichment.csv")

#ZF from combo dataset with GG annotations

tg_genes_signames <- both_res %>% filter(tg_BSnIPRE.class == "pos" | tg_SnIPRE.class == "pos") %>% pull(entrezgene.y)
tg_genes_signames <- tg_genes_signames[!is.na(tg_genes_signames)]
all_genes_tg <- both_res %>% filter(!is.na(entrezgene.y),!is.na(tg_SnIPRE.class)) %>% pull(entrezgene.y)

#Set univeres of genes to those in results for all tests
tg_genes_k <- enrichKEGG(tg_genes_signames,organism="gga",pvalueCutoff=1,pAdjustMethod="BH",qvalueCutoff=1,universe=as.character(all_genes_tg),keyType="ncbi-geneid")

cnetplot(tg_genes_k)
dotplot(tg_genes_k)
summary(tg_genes_k)

write_csv(as.data.frame(tg_genes_k),"TestSpecies_2019/Tguttata_pathwayEnrichment.csv")



#ZF from ZF data alone
gv_genes_signames <- gv_combo_res %>%
  left_join(zf_gene_ids,by=c("gene_name" = "Symbol")) %>%
  filter(BSnIPRE.class == "pos" | SnIPRE.class == "pos") %>%
  filter(!is.na(GeneID)) %>%
  pull(GeneID)
all_genes_gv <- gv_combo_res %>%
  left_join(zf_gene_ids,by=c("gene_name" = "Symbol")) %>%
  filter(!is.na(GeneID),!is.na(SnIPRE.class)) %>%
  pull(GeneID)

#Set univeres of genes to those in results for all tests
gv_genes_k <- enrichKEGG(gv_genes_signames,organism="tgu",pvalueCutoff=1,pAdjustMethod="BH",qvalueCutoff=1,universe=as.character(all_genes_gv),keyType="ncbi-geneid")

cnetplot(gv_genes_k)
dotplot(gv_genes_k)
summary(gv_genes_k)

write_csv(as.data.frame(tg_genes_k),"TestSpecies_2019/Gvarius_pathwayEnrichment.csv")





all_fa_sel_hogs <- fa_res_immfunc %>% filter(BSnIPRE.class == "pos") %>% separate(hog,into = c("HOG","hog")) %>% pull(hog)


all_fa_sel_ncbiIDs <- both_res %>%
  filter(fa_BSnIPRE.class == "pos" | fa_SnIPRE.class == "pos") %>%
  filter(!is.na(entrezgene.y)) %>%
  pull(entrezgene.y)

all_tg_sel_ncbiIDs <- both_res %>%
  filter(tg_BSnIPRE.class == "pos" | tg_SnIPRE.class == "pos", !is.na(all_sel)) %>%
  filter(!is.na(entrezgene.y)) %>%
  pull(entrezgene.y)


all_gv_sel_ncbiIDs <- both_res %>%
  filter(gv_BSnIPRE.class == "pos" | gv_SnIPRE.class == "pos", !is.na(all_sel)) %>%
  filter(!is.na(entrezgene.y)) %>%
  pull(entrezgene.y)


all_genes_signames_tg <- both_res %>%
  filter(all_sel == 1, !is.na(entrezgene.y), !is.na(tg_SnIPRE.class)) %>%
  pull(entrezgene.y)

all_genes_signames_fa <- both_res %>%
  filter(all_sel == 1, !is.na(entrezgene.y), !is.na(fa_SnIPRE.class)) %>%
  pull(entrezgene.y)

all_genes_signames_gv <- both_res %>%
  filter(all_sel == 1, !is.na(entrezgene.y), !is.na(gv_SnIPRE.class)) %>%
  pull(entrezgene.y)

  

perc_overlap_fa <- mean(all_fa_sel_ncbiIDs %in% all_genes_signames_fa)
perc_overlap_tg <- mean(all_tg_sel_ncbiIDs %in% all_genes_signames_tg)
perc_overlap_gv <- mean(all_gv_sel_ncbiIDs %in% all_genes_signames_gv)

rand_sel_fa <- 1:100
for (i in 1:100){
  rand_hogs <- both_res %>% filter(!is.na(entrezgene.y))%>% sample_n(length(all_fa_sel_ncbiIDs)) %>% pull(entrezgene.y)
  rand_sel_fa[i] <- mean(rand_hogs %in% all_genes_signames_fa)
}

rand_sel_tg <- 1:100
for (i in 1:100){
  rand_hogs <- both_res %>% filter(!is.na(entrezgene.y))%>% sample_n(length(all_fa_sel_ncbiIDs)) %>% pull(entrezgene.y)
  rand_sel_tg[i] <- mean(rand_hogs %in% all_genes_signames_tg)
}

rand_sel_gv <- 1:100
for (i in 1:100){
  rand_hogs <- both_res %>% filter(!is.na(entrezgene.y))%>% sample_n(length(all_gv_sel_ncbiIDs)) %>% pull(entrezgene.y)
  rand_sel_gv[i] <- mean(rand_hogs %in% all_genes_signames_gv)
}

ggplot(as.data.frame(rand_sel_gv),aes(x=rand_sel_gv)) + geom_histogram(bins=50) + xlim(c(0,0.7)) +
  geom_segment(data=as.data.frame(perc_overlap_gv),aes(x=perc_overlap_gv, xend = perc_overlap_gv, y = 20, yend =0),col="red",size=2,arrow=arrow(length = unit(0.03, "npc"))) +
  geom_segment(data=as.data.frame(perc_overlap_tg),aes(x=perc_overlap_tg, xend = perc_overlap_tg, y = 20, yend =0),col="red",size=2,arrow=arrow(length = unit(0.03, "npc"))) +
  geom_segment(data=as.data.frame(perc_overlap_fa),aes(x=perc_overlap_fa, xend = perc_overlap_fa, y = 20, yend =0),col="red",size=2,arrow=arrow(length = unit(0.03, "npc"))) +
  xlab("proportion overlap") +
  theme_bw()

ggsave("TestSpecies_2019/gv_tg_fa_species_overlap.png",width=7,height=5)




gga_kegg <- download_KEGG(species = "gga")

paths <- c("gga05164","gga03440","gga03460","gga04620","gga04060","gga04512","gga05168","gga04514","gga03410","gga05132")

pathgenes <- gga_kegg$KEGGPATHID2EXTID %>% tbl_df %>% filter(from %in% paths) %>% rename(entrezgene = to, pathway = from)

comppath_genes <- pathgenes %>% left_join(fa_paml_res,by="entrezgene") %>% filter(!is.na(BSnIPRE.est)) %>% select(pathway,BSnIPRE.est,BSnIPRE.class,FDRPval_busted,FDRPval_m1m2,FDRPval_m8m8a,Omega_m2,Omega_m8,entrezgene)

ggplot(comppath_genes,aes(pathway,BSnIPRE.est)) + geom_violin()

fa_paml_res <- fa_paml_res %>% mutate(pathway="genome") %>% select(pathway,BSnIPRE.est,BSnIPRE.class,FDRPval_busted,FDRPval_m1m2,FDRPval_m8m8a,Omega_m2,Omega_m8,entrezgene) %>% filter(!(entrezgene %in% comppath_genes$entrezgene))

comppath_genome <- fa_paml_res %>% bind_rows(comppath_genes)

ggplot(comppath_genome,aes(pathway,BSnIPRE.est)) + geom_violin()
ggplot(comppath_genome,aes(pathway,BSnIPRE.est)) + geom_boxplot()

comppath_genome %>% filter(!is.na(BSnIPRE.est)) %>% group_by(pathway) %>% summarise(counts = n(), pos_sel = sum(BSnIPRE.class == "pos"), neg_sel = sum(BSnIPRE.class == "neg"), neut = sum(BSnIPRE.class=="neut")) %>% mutate(prop_pos = round(as.numeric(pos_sel)/as.numeric(counts),3)) 



ggplot(comppath_genome,aes(BSnIPRE.est,-log(FDRPval_m1m2),col=pathway)) + geom_point()



for (i in 1:length(paths)){
  print(paths[i])
print(wilcox.test((comppath_genes %>% filter(pathway == paths[i]) %>% pull(BSnIPRE.est)),(fa_res_immfunc %>% pull(BSnIPRE.est))))
}
  


mean(all_fa_sel_hogs %in% all_selected_hogs)


rand_sel <- 1:100
for (i in 1:100){
  rand_hogs <- fa_res_immfunc %>% sample_n(length(all_fa_sel_hogs)) %>% separate(hog,into = c("HOG","hog")) %>% pull(hog)
  rand_sel[i] <- mean(rand_hogs %in% all_selected_hogs)
}



immcols <- c("#e7298a","#7570b3","#1b9e77","#d95f02")
names(immcols) <- c("AllGenes","Receptor","Signaling","Effector")

pdf(file="fa_bayesian_selectioneffect.pdf")
ggplot(data=fa_res_immfunc,aes(Class,BSnIPRE.est)) + geom_violin(aes(fill=Class))
dev.off()







#######Chicken (EB only)
gg_ml_res <- read_csv("gg_eb_results.csv")
gg_rna_to_gene <- read_delim("all_info.summary_galGal",delim="\t",col_names = c("species",	"gene_id","gene_acc","gene_name","biotype","trans_id","trans_acc","trans_len","cds_len","prot_id","is_longest_prot"))
gg_rna_to_gene <- gg_rna_to_gene %>% mutate(gene_acc = as.character(gene_acc))

gg_combo_res <- gg_ml_res %>% select(geneID,PR,FR,PS,FS,Tsil,Trepl,nout,npop,SnIPRE.class,SnIPRE.est,SnIPRE.lbound,SnIPRE.ubound) %>% rename(trans_id = geneID)

#Get geneIDs, only keep longest transcripts (those used in alignment), and drop duplicate proteins, join to get chicken gene IDs and therefore hogs.
gg_combo_res_renamed <- gg_combo_res %>% left_join(fa_rna_to_gene,by = "trans_id") %>% filter(is_longest_prot=="Y") %>% distinct(gene_acc,.keep_all=TRUE) %>% left_join(hog_id_map,by="gene_acc") %>% left_join(hog_id_gg,by="hog") %>% distinct(hog,.keep_all=TRUE)

gg_res_immfunc <- gg_combo_res_renamed %>% left_join(immune,by="entrezgene") %>% mutate(Class = ifelse(is.na(Class),"NonImmune",Class))

gg_res_counts <- gg_res_immfunc %>% group_by(Class) %>% summarise(counts = n(), pos_sel = sum(SnIPRE.class == "pos"), neg_sel = sum(SnIPRE.class == "neg"), neut = sum(SnIPRE.class=="neut")) %>% mutate(prop_pos = round(as.numeric(pos_sel)/as.numeric(counts),3))
sum(gg_res_counts %>% pull(counts))

write_csv(gg_combo_res,path = "gg_mk_results_combined.csv")


all_genes_signames <- gg_res_immfunc %>% filter(SnIPRE.class == "pos") %>% pull(entrezgene)
all_genes_signames <- all_genes_signames[!is.na(all_genes_signames)]
all_genes <- gg_res_immfunc %>% filter(!is.na(entrezgene)) %>% pull(entrezgene)

#Set univeres of genes to those in results for all tests
all_genes_k <- enrichKEGG(all_genes_signames,organism="gga",pvalueCutoff=1,pAdjustMethod="BH",qvalueCutoff=1,universe=as.character(all_genes),keyType="ncbi-geneid")

summary(all_genes_k)

all_gg_sel_hogs <- gg_res_immfunc %>% filter(SnIPRE.class == "pos") %>% separate(hog,into = c("HOG","hog")) %>% pull(hog)

all_gg_sel_ncbiIDs <- gg_res_immfunc %>% filter(SnIPRE.class == "pos") %>% pull(entrezgene)

mean(all_gg_sel_hogs %in% all_selected_hogs)
