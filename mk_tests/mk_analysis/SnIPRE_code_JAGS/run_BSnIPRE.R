#!/usr/bin/Rscript
####
# Run SNIPRE for a given dataset
####

#Example input file:


args = commandArgs(trailingOnly = TRUE)
InputFile = args[1]


#################################################################
## load datasets
#################################################################
source("/n/holylfs/LABS/informatics/ashultz/AvianImmune/popgen/SnIPRE_code_JAGS/B_SnIPRE_source.R")
source("/n/holylfs/LABS/informatics/ashultz/AvianImmune/popgen/SnIPRE_code_JAGS/SnIPRE_source.R")
source("/n/holylfs/LABS/informatics/ashultz/AvianImmune/popgen/SnIPRE_code_JAGS/my.jags2.R")
library(lme4)
library(R2jags)
library(arm)
library(MASS)

data <- read.table(InputFile, header = TRUE)  # sample data set
data$Trepl =  as.integer(data$Trepl)
data$Tsil =  as.integer(data$Tsil)
data <- data[(data$Trepl/data$Tsil)<5,]


#################################################################
## Part (2)  Bayesian Implementation (JAGS package,
##          B_SnIPRE_source.R) necessary
#################################################################




#BSnIPRE.run <- function(mydata, path = ".", burnin = 500, thin = 5, iter = 2500){
# path will be where the chains are stored, and must also be where the ".bug" model is located
# burnin, thin, and iter (number iterations after burnin) are for MCMC samples
BSnIPRE.run(data, burnin = 10000, thin = 4, iter = 15000)

# check to make sure it finished correctly:
# if a "sample" file is in your working directory (getwd()), or the path you sepecified)
# is empty or not there, there is a problem

load("samples")
res.mcmc <- samples
b.res <- BSnIPRE(res.mcmc, data)
bres = b.res$new.dataset
write.table(bres, file = paste(InputFile, ".bayesianresults",sep=""), sep  = ",", row.names = FALSE)