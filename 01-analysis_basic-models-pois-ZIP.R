############################
# Supplemental materials for:
# B. W. Rolek, D. J. Harrison, D. W. Linden,  C. S. Loftin, 
# P. B. Wood. 2021. Associations among breeding 
# conifer-associated birds, forestry treatments, 
# years-since-harvest, and vegetation characteristics in 
# regenerating stands. 
#############################
#############################################
# This file is not needed to replicate analyses.
# We provide the simplest version of the model used
# as a template for other researchers to use for their 
# own analyses
#############################################
# software used
# JAGS 4.3.0 
# R version 4.0.2 (2020-06-22) -- "Taking Off Again"
library (jagsUI) # v1.5.1
load ("./DATA.Rdata")
# set up data
datalfoc$SPP <- length(spp.list.foc)
yr <- array(NA, dim=c(dim (ab)[1], 9) )
yr[,1:3] <- 1; yr[,4:6] <- 2; yr[,7:9] <- 3
datalfoc$yr <- yr
s.year <- array(NA, dim=c(114, 9))
s.year[,1:3] <- 1; s.year[,4:6] <- 2; s.year[,7:9] <- 3
datalfoc$s.year <- s.year

datalfoc$ba <- datalfoc$CovsLam[, "ba"]
nobs <- datalfoc$nobs
dclass <- datalfoc$dclass
int <- datalfoc$int
site <- datalfoc$site
yr_rot <- datalfoc$yr_rot

# print sample sizes 
apply(ab2[,1:2,,,dimnames(ab2)[[5]] %in% spp.list.foc], c(5), sum, na.rm=T)

####################################
# (1) basic Poisson model
####################################
cat("
    model {
    ##### PRIORS ###############################################
    pa.beta <- logit(p.pa.beta0)
    p.pa.beta0 ~ dunif(0,1)
    pp.beta ~ dunif(0, 250)
    lam.beta ~ dnorm(0, 0.01)
    
    ##### DISTANCE AND REMOVAL #####################################
    for (l in 1:L) {
    int[l] ~ dcat(pi.pa.c[site[l], yr_rot[l], ]) # removal class frequencies
    dclass[l] ~ dcat(pi.pd.c[site[l], yr_rot[l], ]) # distance class frequencies
    } # L
    
    # Distance
    for(b in 1:nD){
      f[b] <- (2*midpt[b]*delta)/(B*B)  # radial density function for point counts, change for line transects
    }
    for (i in 1:nsites){
    for(t in 1:YR){  
    for(b in 1:nD){
    g[i,t,b] <- exp(-midpt[b]*midpt[b]/(2*dist.sigma[i,t]*dist.sigma[i,t])) # half-normal distance function
    pi.pd[i,t,b] <- g[i,t,b]*f[b]
    pi.pd.c[i,t,b] <- pi.pd[i,t,b]/pdet[i,t]
    } #nD
    pdet[i,t] <- sum(pi.pd[i,t,1:nD]) # Distance class probabilities
    
    # Removal 
    for (r in 1:R){
    pi.pa[i,t,r] <- p.a[i,t]*pow(1-p.a[i,t], (r-1))
    pi.pa.c[i,t,r] <- pi.pa[i,t,r] / pcap[i,t]
    }  #R
    pcap[i,t] <- sum(pi.pa[i,t,1:R])
    
    # Detection models 
    pmarg[i,t] <-  pcap[i,t]  * pdet[i,t]
    logit(p.a[i,t]) <- pa.beta # add covariates for availability (time-removal) here
    log(dist.sigma[i,t]) <- log(pp.beta) # add covariates for perceptibility (distance) here
    
    ##### POINT-LEVEL ABUNDANCE ###########################     
    nobs[i,t] ~ dbin(pmarg[i,t], N[i,t])  
    N[i,t] ~ dpois(lambda[i,t])
    log(lambda[i,t]) <- lam.beta # add covariates for abundance here
    
    ##### GOODNESS OF FIT #######################################
    nobs.fit[i,t] ~ dbin(pmarg[i,t], N[i,t]) # create new realization of model
    e.p[i,t] <- pmarg[i,t] * N[i,t] # original model prediction
    E.p[i,t] <- pow((nobs[i,t]- e.p[i,t]),2)/(e.p[i,t]+0.5)
    E.New.p[i,t]<- pow((nobs.fit[i,t]-e.p[i,t]),2)/(e.p[i,t]+0.5)
    }} #YR #nsites 
    fit.p <- sum(E.p[1:nsites,1:YR])
    fit.new.p <- sum(E.New.p[1:nsites,1:YR])
    bayesp<-step(fit.new.p-fit.p) # Bayesian p-value for availability model. =0.5 is good fit, near 0 or 1 is poor fit
    
    ##### DERIVED QUANTITIES ####################################
    for(t in 1:YR){
    Ntot[t] <- sum(N[1:nsites,t])
    D[t] <- Ntot[t] / ((3.14*B*B*nsites)/10000)  # dens per ha
    } #YR
    } # End model
    ",file="./model-basic-Poisson.txt")

for (i in 1:19){ 
  spp <- spp.list.foc[i]
  spp.num<- which(dimnames(nobs)[[3]]==spp)
  datalfoc$nobs <- Nav <- apply(ab2[,1:2,,,spp], c(1,4),sum, na.rm=T)
  Mst <- apply(Nav, c(1), max, na.rm=T) +1
  
  inits <- function(){  list(
    N = Nav,
    p.pa.beta0= runif(1, 0.01, 0.99),
    pp.beta= runif(1, 5, 100)
  )  }
  
  params <- c("pa.beta", "p.pa.beta", "pp.beta", 
              "lam.beta",  
              "Ntot", "D",  "bayesp"
  )
  
  # MCMC settings
  ni <- 200000  ;   nb <- 100000   ;   nt <- 10   ;   nc <- 6 ; na=10000
  # Run JAGS
  out <- jags(datalfoc, inits=inits, params, 
              "./model-basic-Poisson.txt", 
              n.thin=nt, n.chains=nc, 
              n.burnin=nb, n.iter=ni, n.adapt=na,
              parallel = T, modules=c("glm")
  )
  
  fn<- paste( "./", spp, "_basic-pois.RData", sep="" )
  save(list= c("out"), file=fn)
}

###################################
# (2) basic zero-inflated Poisson model
###################################
cat("
    model {
    ##### PRIORS ###############################################
    pa.beta <- logit(p.pa.beta)
    p.pa.beta ~ dunif(0,1)
    pp.beta ~ dunif(0, 250)
    lam.beta ~ dnorm(0, 0.01)
    psi.beta <- logit(p.psi.beta)
    p.psi.beta ~ dunif(0,1)
    
    ##### DISTANCE AND REMOVAL #####################################
    for (l in 1:L) {
    int[l] ~ dcat(pi.pa.c[site[l], yr_rot[l], ]) # removal class frequencies
    dclass[l] ~ dcat(pi.pd.c[site[l], yr_rot[l], ]) # distance class frequencies
    } # L
    
    # Distance 
    for(b in 1:nD){
      f[b] <- (2*midpt[b]*delta)/(B*B)  # radial density function for point counts, change for line transects
    }
    for (i in 1:nsites){
    for(t in 1:YR){  
    for(b in 1:nD){
    g[i,t,b] <- exp(-midpt[b]*midpt[b]/(2*dist.sigma[i,t]*dist.sigma[i,t])) # half-normal distance function
    pi.pd[i,t,b] <- g[i,t,b]*f[b]
    pi.pd.c[i,t,b] <- pi.pd[i,t,b]/pdet[i,t]
    } #nD
    pdet[i,t] <- sum(pi.pd[i,t,1:nD]) # Distance class probabilities
    
    # Removal 
    for (r in 1:R){
    pi.pa[i,t,r] <- p.a[i,t]*pow(1-p.a[i,t], (r-1))
    pi.pa.c[i,t,r] <- pi.pa[i,t,r] / pcap[i,t]
    }  #R
    pcap[i,t] <- sum(pi.pa[i,t,1:R])
    
    # Detection models 
    pmarg[i,t] <-  pcap[i,t]  * pdet[i,t]
    logit(p.a[i,t]) <- pa.beta # add covariates for availability (time-removal) here
    log(dist.sigma[i,t]) <- log(pp.beta) # add covariates for perceptibility (distance) here
    
    ##### POINT-LEVEL ABUNDANCE ###########################     
    nobs[i,t] ~ dbin(pmarg[i,t], N[i,t])  
    N[i,t] ~ dpois(lambda.eff[i,t])
    lambda.eff[i,t] <- lambda[i,t] * w.lam[i,t]
    w.lam[i,t] ~  dbern(psi[i,t])
    log(lambda[i,t]) <- lam.beta # add covariates for abundance
    logit(psi[i,t]) <- psi.beta # add covariates for habitat suitability
    # If not running set inits near psi=1 on logit scale. For example psi.beta=10
    
    ##### GOODNESS OF FIT #######################################
    nobs.fit[i,t] ~ dbin(pmarg[i,t], N[i,t]) # create new realization of model
    e.p[i,t] <- pmarg[i,t] * N[i,t] # original model prediction
    E.p[i,t] <- pow((nobs[i,t]- e.p[i,t]),2)/(e.p[i,t]+0.5)
    E.New.p[i,t]<- pow((nobs.fit[i,t]-e.p[i,t]),2)/(e.p[i,t]+0.5)
    }} #YR #nsites 
    fit.p <- sum(E.p[1:nsites,1:YR])
    fit.new.p <- sum(E.New.p[1:nsites,1:YR])
    bayesp<-step(fit.new.p-fit.p) # Bayesian p-value for availability model. =0.5 is good fit, near 0 or 1 is poor fit
    
    ##### DERIVED QUANTITIES ####################################
    for(t in 1:YR){
    Ntot[t] <- sum(N[1:nsites,t])
    D[t] <- Ntot[t] / ((3.14*B*B*nsites)/10000)  # dens per ha
    } #YR
    } # End model
    ",file="./model-basic-ZIP.txt")

for (i in 1:19){ 
  spp <- spp.list.foc[i]
  spp.num<- which(dimnames(nobs)[[3]]==spp)
  datalfoc$nobs <- Nav <- apply(ab2[,1:2,,,spp], c(1,4),sum, na.rm=T)
  Mst <- apply(Nav, c(1), max, na.rm=T) +1
  
  inits <- function(){  list(
    N = Nav, p.psi.beta= 0.99, # setting psi near 1 helps model run
    p.pa.beta0= runif(1, 0.1, 0.9),
    pp.beta= runif(1, 5, 100)
  )  }
  
  params <- c("pa.beta", "p.pa.beta", "pp.beta", 
              "lam.beta", "psi.beta", "p.psi.beta",
              "Ntot", "D",  "bayesp"
  )
  
  # MCMC settings
  ni <- 200000  ;   nb <- 100000   ;   nt <- 10   ;   nc <- 6 ; na=10000
  # Run JAGS
  out <- jags(datalfoc, inits=inits, params, 
               "./model-basic-ZIP.txt", 
              n.thin=nt, n.chains=nc, 
              n.burnin=nb, n.iter=ni, n.adapt=na,
              parallel = T, modules=c("glm")
              )
  
  fn<- paste( "./", spp, "_basic-ZIP.RData", sep="" )
  save(list= c("out"), file=fn)
}

