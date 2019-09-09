# G(n) from equation 1.4 in dissertation

Gnint <- function(s,n=30){
int <- 0
  f <- function(x,g){(1-x)^(n-1) * ( 1- exp(-2*g*x))/(2*g*x)}
for(i in 1:length(s)){
  int[i] <- integrate(f,g= s[i], 0,1)$value}
return(int)
}

# F(n) from equation 1.3 in dissertation

Fnint <- function(s, n = 30){
int <- 0
  f <- function(x,g){(1-x^(n) - (1-x)^(n))/(1-x)* ( 1- exp(-2*g*x))/(2*g*x)}
for(i in 1:length(s)){
  int[i] <- integrate(f,g= s[i], 0,1)$value}
return(int)
} 

# solution to "func" is the selection effect 

func <- function(g, ef, tau, n=30, m =1){
  L.n <- sum(1/(c(1:(n-1))))
  f <- ef - log((tau +Gnint(g,m)+Gnint(g,n))/(Fnint(g, n))*((L.n)/(tau+1/m+1/n)))
return(f)
}

# uses the unit root function to find root of "func" (the estimated selection effect)
    
LLest2gamm <- function(effect,tau.est = rep(10, length(effect)), 
                       n = rep(30, length(effect)),m =rep(1, length(effect))){
gest <- 0
neg <- which(tau.est < 0)
tau.est[neg] <- 0
for(i in 1:length(effect)){
    worked <- try(gest[i] <- uniroot(func, lower=-350, upper=300,
                                     ef=effect[i],tau=tau.est[i], n=n[i], m = m[i])$root)
    if (class(worked)=="try-error") gest[i] <- NA
                }
return(gest)
}


# estimate mutation rate

theta <- function(gene.effect.est, m, n){
  th <- exp(gene.effect.est)/(Ln(m)+Ln(n))
  return(th)
}

# function calls estimation of mutation rate, and estimates tau 
tau.theta <- function(gene.effect.est,f.effect.est,
                      m=rep(1, length(gene.effect.est)), 
                      n=rep(30, length(gene.effect.est))){
  th <- theta(gene.effect.est, m, n)
  t <- exp(gene.effect.est+f.effect.est)/th  - 1/m - 1/n
  return(list(tau.est = t, theta.est=th))
}


# root of function is estimate f (g = estimate of gamma)
funcf <- function(frac, g, r.eff, npop=30, nout=1){
  denom <- Ln(npop)+Ln(nout)
  f <- r.eff - log(frac*(2*g)/(1-exp(-2*g))*(Fnint(g,npop)+Fnint(g, nout))/denom)
  return(f)
}

#functions that find roots of funcf
frac.est <- function(r.eff, g, npop = rep(30, length(r.eff)), 
                    nout = rep(1, length(r.eff))){
  frac <- 0
  for(i in 1:length(r.eff)){

    worked <- try(frac[i] <-  uniroot(funcf, lower = 0, upper = 100, 
                       r.eff = r.eff[i],g = g[i], npop=npop[i], nout=nout[i])$root)
    if (class(worked)=="try-error") frac[i] <- NA

  }
  return(frac)
}

# L(n) from equation 1.2 in disseration
Ln <- function(n){
  v <- 0
  t <- 0
  for(i in 1:length(n)){
  v <- c(1:(n[i]-1))
  t[i] <- sum(1/v)
  if(n[i]==1) t[i] <- 0}
  return(t)
}

# calculates the probability of being "negative"

pn <- function(vec){
  negprop <- length(which(vec < 0))/length(vec)
  return(negprop)}
# prep and run WinBUGS
BSnIPRE.run <- function(mydata, path = ".", burnin = 500, thin = 5, iter = 2500){
  #path will be where the chains are stored, and must also be where the ".bug" model is located
data <- mydata
PS <- data$PS
PR <- data$PR
FR <- data$FR
FS <- data$FS

# subset of data
n <- length(FS) #number of genes
N <- 4*n #number of observations
TS <- data$Tsil
TR <- data$Trepl

Ivec <- matrix(1, nrow = n)  # makes one vector of appropriate size subset
d.mu <-as.numeric( matrix(c(PS,PR,FS,FR), ncol =1))
d.replacement <- as.numeric(c(0,1,0,1)%x%Ivec)
d.fixed <- as.numeric(c(0,0,1,1)%x%Ivec)
d.gene <- as.vector(rep(1,4)%x%c(1:n))
d.TS <- as.vector(rep(1,4)%x%c(TS))
d.TR <- as.vector(rep(1,4)%x%c(TR))

m = c(0,0,0,0)
count <- as.vector(d.mu)
R <- as.vector(d.replacement)
F <- as.vector(d.fixed)
G <- as.vector(d.gene)
RF <- as.vector(R*F)
x <- matrix(c(rep(1,4*n),R,F,RF), nrow = 4*n)
TR <- as.vector(d.TR)
TS <- as.vector(d.TS)
v = 10  # df 

iS <- diag(4)
S <- solve(iS)/v

wb.sfs1 <- list("count"=count,"R"=R,"n"=n,"N"=N,"G"=G, "TR"=TR, "TS"=TS,"x"=x,"S"=S)
inits = function(){list(
  BG = array(mvrnorm(n,c(0,0,0,0),iS),c(n,4)),
  T = solve(iS),
  m = c(0,0,0,0)
 )}

my.inits = list(inits(), inits(),inits())

setwd(path)

constr.sim <- my.jags(data=wb.sfs1,inits=my.inits,  n.chains = 3,
                    model.file = "SnIPRE.bug",      # v = 10
                     parameters.to.save = c("m","BG","T"),
#                   program = "OpenBUGS",
                   n.iter = iter, n.thin=thin, n.burnin = burnin)

}

## Main function call for B SnIPRE results analysis
BSnIPRE <- function(data.mcmc,mydata){
  n = dim(mydata)[1]
  npop = mydata$npop
  nout = mydata$nout
  BRF <- data.mcmc[,which(varnames(data.mcmc) == "m[4]")]
  BR <- data.mcmc[,which(varnames(data.mcmc) == "m[2]")]
  BF <- data.mcmc[,which(varnames(data.mcmc)=="m[3]")]
  B <- data.mcmc[,which(varnames(data.mcmc)=="m[1]")]
  BGst <- which(varnames(data.mcmc) == "BG[1,1]")
  BG <-data.mcmc[, c(BGst:(BGst+4*n-1))]
  
  RF <- c(3*n+(1:n))  #4
  R <- c(n+(1:n))        #2
  F <- c(2*n+(1:n))        #3
  G <- c(1:n)        #1

  gene.effect1 <- BG[,G,][[1]]
  gene.effect2 <- BG[,G,][[2]]
  gene.effect3 <- BG[,G,][[3]]

  gene.effect <- mcmc.list(gene.effect1, gene.effect2, gene.effect3)
  
  sel.effect1 <- BG[,RF,][[1]]
  sel.effect2 <- BG[,RF,][[2]]
  sel.effect3 <- BG[,RF,][[3]]

  sel.effect <- mcmc.list(sel.effect1, sel.effect2, sel.effect3)

  r.effect1 <- BG[,R,][[1]]
  r.effect2 <- BG[,R,][[2]]
  r.effect3 <- BG[,R,][[3]]

  r.effect <- mcmc.list(r.effect1, r.effect2, r.effect3)

  f.effect1 <- BG[,F,][[1]]
  f.effect2 <- BG[,F,][[2]]
  f.effect3 <- BG[,F,][[3]]

  f.effect <- mcmc.list(f.effect1, f.effect2, f.effect3)
  
 
  sel.probs <- apply(as.matrix(sel.effect), 2, quantile, probs = c(.025, .5, .975))
  pneg <- apply(as.matrix(sel.effect), 2, pn)
  r.probs <- apply(as.matrix(r.effect), 2, quantile, probs = c(.025, .5, .975))
                                        # Selection Effect Estimate/Classification
  mydata$BSnIPRE.class <- 0
  mydata$BSnIPRE.est <- sel.probs[2,]
  mydata$BSnIPRE.lbound <- sel.probs[1,]
  mydata$BSnIPRE.ubound <- sel.probs[3,]
  mydata$BSnIPRE.class = "neut"
  mydata$BSnIPRE.class[which(mydata$BSnIPRE.lbound>0)] = "pos"
  mydata$BSnIPRE.class[which(mydata$BSnIPRE.ubound<0)] = "neg"
  mydata$BSnIPRE.pneg <- pneg

                                        # Replacement Effect Estimate/Classification
  mydata$BSnIPRE.Rest <- r.probs[2,]
  mydata$BSnIPRE.Rlbound <- r.probs[1,]
  mydata$BSnIPRE.Rubound <- r.probs[3,]
  
  mydata$BSnIPRE.Rclass <- 0
  mydata$BSnIPRE.Rclass = "neut"
  mydata$BSnIPRE.Rclass[which(mydata$BSnIPRE.Rlbound>0)] = "pos"
  mydata$BSnIPRE.Rclass[which(mydata$BSnIPRE.Rubound<0)] = "neg"

                                        # gene effect estimate
 probs.gene <- apply(as.matrix(gene.effect), 2, quantile, 
                     probs = c(.025,.5,.975))
                                        # fixed effect estimate 
 probs.fix <- apply(as.matrix(f.effect), 2, quantile, probs = c(.025,.5,.975))   
 gene.effect.est <- probs.gene[2,]
 f.effect.est <- probs.fix[2,]
                                        # estimate tau & theta
 params <- tau.theta(gene.effect.est,f.effect.est, npop, nout)     
 mydata$BSnIPRE.tau = params$tau.est
 mydata$BSnIPRE.theta = params$theta.est
 mydata$BSnIPRE.gamma[mydata$BSnIPRE.est>-4.35] <-  LLest2gamm(
                   mydata$BSnIPRE.est[mydata$BSnIPRE.est>-4.35],
                   mydata$BSnIPRE.tau[mydata$BSnIPRE.est>-4.35],
                   n = npop, m = nout)
na.set <- which(is.na(mydata$BSnIPRE.gamma) == TRUE)
good.set <- which(mydata$BSnIPRE.est > - 4.4)
use <- setdiff(good.set, na.set)
mydata$BSnIPRE.f <- NA
mydata$BSnIPRE.f.lb <- NA
mydata$BSnIPRE.f.ub <- NA
g.ests <- mydata$BSnIPRE.gamma
g.zeros <- which(mydata$BSnIPRE.class == "neut")
g.ests[g.zeros] <- .000001
mydata$BSnIPRE.f[use] <- frac.est(mydata$BSnIPRE.Rest[use], g.ests[use],
                                          npop[use], nout[use])
mydata$BSnIPRE.f.lb[use] <-frac.est(mydata$BSnIPRE.Rlbound[use],
                                   g.ests[use],  npop[use],nout[use])
mydata$BSnIPRE.f.ub[use] <-frac.est(mydata$BSnIPRE.Rubound[use], 
                                   g.ests[use], npop[use],nout[use])
  mydata$BSnIPRE.f.class <- "neut"
  mydata$BSnIPRE.f.class[which(mydata$BSnIPRE.f.ub<1)] = "neg"
  mydata$BSnIPRE.f.class[which(mydata$BSnIPRE.f.lb>1)] = "pos"
return(list(new.dataset = mydata, sel.effect = sel.effect, rep.effect = r.effect,
              fix.effect = f.effect, gene.effect = gene.effect) )
}
