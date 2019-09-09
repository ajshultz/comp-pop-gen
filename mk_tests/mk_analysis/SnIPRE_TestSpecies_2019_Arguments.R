#!/usr/bin/env Rscript

## Part (1)  Empirical Bayes Implementation  (lme4 package, SnIPRE_source.R)
## Part (2)  Bayesian Implementation (R2WinBUGS package, B_SnIPRE_source.R, and WinBUGS or OpenBUGS) necessary
setwd("~/Dropbox/BirdImmuneGeneEvolution/PopGen/TestSpecies_2019")

#################################################################
## Part (1)  Empirical Bayes Implementation  (lme4 package)
#################################################################

source("SnIPRE_code_JAGS/SnIPRE_source.R")
source("SnIPRE_code_JAGS/my.jags2.R")
library(lme4)
library(R2jags)
library(arm)
library(tidyverse)

#Arguments for analyses, first should be species identification, second should be the number of indiviuals squenced and third should be the number of outgroup individuals

args = commandArgs(trailingOnly = TRUE)
species = args[1]
npop_sp = args[2]
nout_sp = args[3]

#species="Tguttata"

#species="Gvarius"
#nout_sp=2
#npop_sp=9

#Read Tguttata data
tg_data <- read_delim(paste0(species,"_PolymorphismDivergenceStats_combined.txt"),delim="\t")
#Format to fit SnIPRE requirements (geneID, PR (# nonsynpoly), FR (#nonsyndiv), PS (# synpoly), FS (# syndiv), Tsil (# syn sites), Trepl (# nonsyn sites), nout (# outgroup seqs), npop (# ingroup seqs))
#For now just putting the number of inds sequenced as # outgroup seqs (2 species, 21 inds), and # ingroup seqs (19 inds)
tg_data <- tg_data %>%
  dplyr::select(geneID = Gene, PR = NonsynPolymorphic, FR = NonsynFixed, PS = SynPolymorphic, FS = SynFixed, Tsil = CallableSynSites, Trepl = CallableNsynSites) %>%
  mutate(nout = nout_sp, npop = npop_sp) %>%
  mutate_at(vars(PR:npop),.fun=as.integer)

tg_data <- tg_data %>%
  filter((Trepl/Tsil)<5) %>%
  filter((PR+FR+PS+FS)>1) %>%
  as.data.frame

#SnIPRE <-function(mydata)
# mydata: name of data set;
# mydata must have a header with the following columns: PS, PR, FS, FR, npop, nout, Tsil, Trepl (no particular order)
# outputs 2 objects:  new.dataset & model
# new.dataset contains the estimates of the selection effect, and selection coefficient (gamma); as well as the estimates of constraint effect (Rest) and constraint (f). 
tg.res = SnIPRE(tg_data)

tg.qres = tg.res$new.dataset
tg.model = tg.res$model

table(tg.qres$SnIPRE.class)

write.table(tg.qres, file = paste0(species,"_results.csv"), sep  = ",", row.names = FALSE)






#################################################################
## Part (2)  Bayesian Implementation (R2WinBUGS package,
##          B_SnIPRE_source.R, and WinBUGS or OpenBUGS) necessary
#################################################################
setwd("SnIPRE_code_JAGS")
source("B_SnIPRE_source.R")
source("my.jags2.R")
library(lme4)
library(R2jags)
library(arm)


#data <- read.table("data_example.txt", header = TRUE)  # sample data set
#BSnIPRE.run <- function(mydata, path = ".", burnin = 500, thin = 5, iter = 2500){
  # path will be where the chains are stored, and must also be where the ".bug" model is located
  # burnin, thin, and iter (number iterations after burnin) are for MCMC samples


BSnIPRE.run(tg_data, burnin = 500, thin = 4, iter = 2500)

# check to make sure it finished correctly:
# if a "sample" file is in your working directory (getwd()), or the path you sepecified)
# is empty or not there, there is a problem


load("samples")

res.mcmc <- samples

#BSnIPRE <- function(data.mcmc,mydata){
# outputs 2 objects:  new.dataset & effects
# new.dataset contains the estimates of the selection effect, and selection coefficient (gamma); as well as the estimates of constraint effect (Rest) and constraint (f). 
# the "effects" may be useful if you are interested in estimation
# of population parameters (gamma, constraint) with other assumptions than the PRF

tg.b.res <- BSnIPRE(res.mcmc, tg_data)

tg.bres = tg.b.res$new.dataset

write.table(tg.bres, file = paste0(species,"_bayesian_results.csv"), sep  = ",", row.names = FALSE)



