## Part (1)  Empirical Bayes Implementation  (lme4 package, SnIPRE_source.R)
## Part (2)  Bayesian Implementation (R2WinBUGS package, B_SnIPRE_source.R, and WinBUGS or OpenBUGS) necessary
setwd("~/Dropbox/SnIPRE_code_JAGS")

#################################################################
## Part (1)  Empirical Bayes Implementation  (lme4 package)
#################################################################

source("SnIPRE_source.R")
source("my.jags2.R")
library(lme4)
library(R2jags)
library(arm)



data <- read.table("data_example.txt", header = TRUE)  # sample data set

#SnIPRE <-function(mydata)
# mydata: name of data set;
# mydata must have a header with the following columns: PS, PR, FS, FR, npop, nout, Tsil, Trepl (no particular order)
# outputs 2 objects:  new.dataset & model
# new.dataset contains the estimates of the selection effect, and selection coefficient (gamma); as well as the estimates of constraint effect (Rest) and constraint (f). 
eb.res = SnIPRE(data)

res = eb.res$new.dataset
model = eb.res$model

write.table(res, file = "eb_results.csv", sep  = ",", row.names = FALSE)



#################################################################
## Part (2)  Bayesian Implementation (R2WinBUGS package,
##          B_SnIPRE_source.R, and WinBUGS or OpenBUGS) necessary
#################################################################

source("B_SnIPRE_source.R")
source("my.jags2.R")
library(lme4)
library(R2jags)
library(arm)


data <- read.table("data_example.txt", header = TRUE)  # sample data set
#BSnIPRE.run <- function(mydata, path = ".", burnin = 500, thin = 5, iter = 2500){
  # path will be where the chains are stored, and must also be where the ".bug" model is located
  # burnin, thin, and iter (number iterations after burnin) are for MCMC samples
BSnIPRE.run(data, burnin = 10000, thin = 4, iter = 15000)

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

b.res <- BSnIPRE(res.mcmc, data)

bres = b.res$new.dataset

write.table(bres, file = "bayesian_results.csv", sep  = ",", row.names = FALSE)

