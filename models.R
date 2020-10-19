# This file includes all model files used
# for Rolek et al. 202X
# Includes 6 distance-removal abundance models
# to estimate abundance.
# basic models- Poisson (1) and ZIP (2)
# vegetation models- Poisson (3) and ZIP (4) 
# treatment models- Poisson (5) and ZIP (6) 

##### Notation ##########################################
## indices: i=site, k=visit, t=year
## XX.beta = parameters for regression
## pa... related to availability
## pp... related to perceptibility
## lam... related to abundance
## psi... related to habitat suitability
## dist.sigma = distance scale parameter
## N = detection corrected abundance
## Ntot = population size of total area surveyed
## D = density
## bayesp = Bayesian p-value for model fit
## w... = Bernoulli indicator variables for GVS

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
    for (i in 1:nsites){
    for(t in 1:YR){  
    for(b in 1:nD){
    g[i,t,b] <- exp(-midpt[b]*midpt[b]/(2*dist.sigma[i,t]*dist.sigma[i,t])) # half-normal distance function
    f[b] <- (2*midpt[b]*delta)/(B*B)     # radial density function for point counts, change for line transects
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
    ",file="./basic-model-Poisson.txt")

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
    for (i in 1:nsites){
    for(t in 1:YR){  
    for(b in 1:nD){
    g[i,t,b] <- exp(-midpt[b]*midpt[b]/(2*dist.sigma[i,t]*dist.sigma[i,t])) # half-normal distance function
    f[b] <- (2*midpt[b]*delta)/(B*B)     # radial density function for point counts, change for line transects
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
    ",file="./basic-model-ZIP.txt")

###############################################
# (3) Vegetation model- 
# Poisson distribution - 
# Gibbs variable selection
###############################################
cat("
    model {
    ##### PRIORS ###############################################
    pa.beta[1] <- logit(p.pa.beta0)
    p.pa.beta0 ~ dunif(0,1) 
    pp.beta[1] ~ dunif(0, 250)
    
    # priors for the w model inclusion terms
    w[13] ~ dbern(0.1667)
    w[12] ~ dbern(0.1667)
    w12_13 <- w[12]+w[13]
    p.w11 <- equals(w12_13,0)*.3077 + (1-equals(w12_13,0))
    w[11] ~ dbern(p.w11)
    p.w10 <- equals(w[13],0)*0.4 + w[13]
    w[10] ~ dbern(p.w10)
    p.w9 <- equals(w[12],0)*0.4 + w[12]
    w[9] ~ dbern(p.w9)
    w[8] ~ dbern(0.5)
    w[7] ~ dbern(0.5)
    w[6] ~ dbern(0.5)
    w[5] ~ dbern(0.5)
    w[4] ~ dbern(0.5)
    w10_13 <- w[10] + w[11] + w[12] + w[13]
    p.w3 <- equals(w10_13,0)*0.5 + (1-equals(w10_13,0))
    w[3] ~ dbern(p.w3)
    w9.11_13 <- w[9] + w[11] + w[12] + w[13]
    p.w2 <- equals(w9.11_13,0)*0.5 + (1-equals(w9.11_13,0))
    w[2] ~ dbern(p.w2)
    w[1] ~ dbern(1)
    
    # priors for wpa for availability
    wpa[6] ~ dbern(0.125)
    p.wpa5 <- (1-wpa[6])*0.143 + wpa[6] 
    wpa[5] ~ dbern(p.wpa5) 
    p.wpa4 <- (1-wpa[6])*0.286 + wpa[6] 
    wpa[4] ~ dbern(p.wpa4)
    wpa456 <- wpa[4]+wpa[5] + wpa[6]
    p.wpa3 <- equals(wpa456,0)*0.5 + (1-equals(wpa456,0)) 
    wpa[3] ~ dbern(p.wpa3)
    wpa56  <- wpa[5]+wpa[6]
    p.wpa2 <- equals(wpa56,0)*0.5 + (1-equals(wpa56,0)) 
    wpa[2] ~ dbern(p.wpa2)
    wpa[1] ~ dbern(1)
    
    # priors for wpp for perceptibility
    wpp[1] ~ dbern(1) #intercept
    for (n in 2:4){wpp[n] ~ dbern(0.5)  }
    
    # set up the vectors/matrices for beta estimation, abundance
    for(b1 in 1:n.betas){
    wtemp[b1] <- w[pos[b1]]   # this uses GVS
    #wtemp[b1] <- 1   # this forces you to fit the full model (all terms included)
    mean.b[b1] <- post.b[b1]*(1-wtemp[b1])  # prior is either 0 or full-model posterior mean
    for(b2 in 1:n.betas){   # set up the precision matrix (inverse variance) # allows for betas to be multivariate, if desired
    tau.b[b1,b2] <- equals(b1,b2)*((1/sd.b[b1]^2)*(1-wtemp[b1])) + (wtemp[b1]*b.tau)
    } # b2
    lam.beta[b1] ~ dnorm(mean.b[b1],tau.b[b1,b1])   # all beta coefficients
    } # b1
    
    # set up the vectors/matrices for beta estimation, availability
    for(b1 in 2:n.betas.pa){ # starts at 2 because intercept always requires a diff prior and w=1
    wpa.temp[b1] <- wpa[pos.pa[b1]]
    #wpa.temp[b1] <- 1
    mean.b.pa[b1] <- post.b.pa[b1]*(1-wpa.temp[b1])
    for(b2 in 2:n.betas.pa){ # starts at 2 because intercept always requires a diff prior and w=1
    tau.b.pa[b1,b2] <- equals(b1,b2)*((1/sd.b.pa[b1]^2)*(1-wpa.temp[b1])) + (wpa.temp[b1]*b.tau.pa)
    } # b2
    pa.beta[b1] ~ dnorm(mean.b.pa[b1],tau.b.pa[b1,b1])   # all beta coefficients
    } # b1
    
    # set up the vectors/matrices for beta estimation, perceptility
    for(b1 in 2:n.betas.pp){ # starts at 2 because intercept always requires a diff prior and w=1
    wpp.temp[b1] <- wpp[pos.pp[b1]]
    #wpp.temp[b1] <- 1
    mean.b.pp[b1] <- post.b.pp[b1]*(1-wpp.temp[b1])
    for(b2 in 2:n.betas.pp){ # starts at 2 because intercept always requires a diff prior and w=1
    tau.b.pp[b1,b2] <- equals(b1,b2)*((1/sd.b.pp[b1]^2)*(1-wpp.temp[b1])) + (wpp.temp[b1]*b.tau.pp)
    } # b2
    pp.beta[b1] ~ dnorm(mean.b.pp[b1],tau.b.pp[b1,b1])   # all beta coefficients
    } # b1
    
    # vector of stand-specific predictors
    lam.mu <- mm[,] %*% (lam.beta[]*wtemp[])
    
    stand.tau <- 1/ (stand.sig*stand.sig)
    stand.sig ~ dunif(0,10)
    yr.tau <- 1/ (yr.sig*yr.sig)
    yr.sig ~ dunif(0,20)
    obs.tau <- 1/ (obs.sig*obs.sig)
    obs.sig ~ dunif(0,20)
    b.tau <- 0.01
    b.tau.pa <- 0.01
    b.tau.pp <- 0.01
    
    ##### DISTANCE AND REMOVAL #####################################
    for (l in 1:L) {
    int[l] ~ dcat(pi.pa.c[site[l], yr_rot[l], ]) # removal class frequencies
    dclass[l] ~ dcat(pi.pd.c[site[l], yr_rot[l], ]) # distance class frequencies
    } # L
    
    # Distance 
    for (i in 1:nsites){
    for(t in 1:YR){  
    for(b in 1:nD){
    g[i,t,b] <- exp(-midpt[b]*midpt[b]/(2*dist.sigma[i,t]*dist.sigma[i,t])) # half-normal distance function
    f[b] <- (2*midpt[b]*delta)/(B*B)     # radial density function for point counts, change for line transects
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
    logit(p.a[i,t]) <- wpa[1]*pa.beta[1] + wpa.temp[2]*pa.beta[2]*hr[i,t] + wpa.temp[3]*pa.beta[3]*date[i,t] + 
    wpa.temp[4]*pa.beta[4]*date2[i,t] +
    wpa.temp[5]*pa.beta[5]*date[i,t]*hr[i,t] + wpa.temp[6]*pa.beta[6]*date[i,t]*date2[i,t]*hr[i,t] 
    log(dist.sigma[i,t]) <- wpp[1]*log(pp.beta[1]) +  wpp.temp[2]*pp.beta[2]*densiom[i,t] + 
    wpp.temp[3]*pp.beta[3]*noise[i,t] + wpp.temp[4]*pp.beta[4]*ba[i] + obs.eps[obs[i,t]]
    
    ##### POINT-LEVEL ABUNDANCE ###########################     
    nobs[i,t] ~ dbin(pmarg[i,t], N[i,t])  
    log(lambda[i,t]) <- lam.beta.s[stand.id[i]] + yr.eps[t] + lam.mu[i]
    N[i,t] ~ dpois(lambda[i,t])
    
    ##### GOODNESS OF FIT #######################################
    nobs.fit[i,t] ~ dbin(pmarg[i,t], N[i,t]) # create new realization of model
    e.p[i,t] <- pmarg[i,t] * N[i,t] # original model prediction
    E.p[i,t] <- pow((nobs[i,t]- e.p[i,t]),2)/(e.p[i,t]+0.5)
    E.New.p[i,t]<- pow((nobs.fit[i,t]-e.p[i,t]),2)/(e.p[i,t]+0.5)
    }} #YR #nsites 
    fit.p <- sum(E.p[1:nsites,1:YR])
    fit.new.p <- sum(E.New.p[1:nsites,1:YR])
    bayesp<-step(fit.new.p-fit.p) # Bayesian p-value for availability model. =0.5 is good fit, near 0 or 1 is poor fit

    # Random effects
    for (s in 1:S) { lam.beta.s[s] ~ dnorm(0, stand.tau) } #S
    for (y in 1:9){ yr.eps[y] ~ dnorm(0, yr.tau)}
    for (o in 1:28){ obs.eps[o] ~ dnorm(0, obs.tau)} 
    
    ##### DERIVED QUANTITIES ####################################
    for(t in 1:YR){
    Ntot[t] <- sum(N[1:nsites,t])
    D[t] <- Ntot[t] / ((3.14*B*B*nsites)/10000)  # dens per ha
    } #YR
    } # End model
    ",file="./veg-Poisson-GVS.txt")

###############################################
# (4) Vegetation model- 
# Zero-inflated Poisson distribution - 
# Gibbs variable selection
###############################################
cat("
    model {
    ##### PRIORS ###############################################
    pa.beta[1] <- logit(p.pa.beta0)
    p.pa.beta0 ~ dunif(0,1) 
    pp.beta[1] ~ dunif(0, 250)
    
    # priors for the w model inclusion terms
    w[13] ~ dbern(0.1667)
    w[12] ~ dbern(0.1667)
    w12_13 <- w[12]+w[13]
    p.w11 <- equals(w12_13,0)*.3077 + (1-equals(w12_13,0))
    w[11] ~ dbern(p.w11)
    p.w10 <- equals(w[13],0)*0.4 + w[13]
    w[10] ~ dbern(p.w10)
    p.w9 <- equals(w[12],0)*0.4 + w[12]
    w[9] ~ dbern(p.w9)
    w[8] ~ dbern(0.5)
    w[7] ~ dbern(0.5)
    w[6] ~ dbern(0.5)
    w[5] ~ dbern(0.5)
    w[4] ~ dbern(0.5)
    w10_13 <- w[10] + w[11] + w[12] + w[13]
    p.w3 <- equals(w10_13,0)*0.5 + (1-equals(w10_13,0))
    w[3] ~ dbern(p.w3)
    w9.11_13 <- w[9] + w[11] + w[12] + w[13]
    p.w2 <- equals(w9.11_13,0)*0.5 + (1-equals(w9.11_13,0))
    w[2] ~ dbern(p.w2)
    w[1] ~ dbern(1)
    
    # priors for wpa for availability
    wpa[6] ~ dbern(0.125)
    p.wpa5 <- (1-wpa[6])*0.143 + wpa[6] 
    wpa[5] ~ dbern(p.wpa5) 
    p.wpa4 <- (1-wpa[6])*0.286 + wpa[6] 
    wpa[4] ~ dbern(p.wpa4)
    wpa456 <- wpa[4]+wpa[5] + wpa[6]
    p.wpa3 <- equals(wpa456,0)*0.5 + (1-equals(wpa456,0)) 
    wpa[3] ~ dbern(p.wpa3)
    wpa56  <- wpa[5]+wpa[6]
    p.wpa2 <- equals(wpa56,0)*0.5 + (1-equals(wpa56,0)) 
    wpa[2] ~ dbern(p.wpa2)
    wpa[1] ~ dbern(1)
    
    # priors for wpp for perceptibility
    wpp[1] ~ dbern(1) #intercept
    for (n in 2:4){wpp[n] ~ dbern(0.5)  }
    
    # set up the vectors/matrices for beta estimation, abundance
    for(b1 in 1:n.betas){
    wtemp[b1] <- w[pos[b1]]                # this uses GVS
    #wtemp[b1] <- 1                          # this forces you to fit the full model (all terms included)
    mean.b[b1] <- post.b[b1]*(1-wtemp[b1])  # prior is either 0 or full-model posterior mean
    for(b2 in 1:n.betas){                   # set up the precision matrix (inverse variance) # allows for betas to be multivariate, if desired
    tau.b[b1,b2] <- equals(b1,b2)*((1/sd.b[b1]^2)*(1-wtemp[b1])) + (wtemp[b1]*b.tau)
    } # b2
    lam.beta[b1] ~ dnorm(mean.b[b1],tau.b[b1,b1])   # all beta coefficients
    } # b1
    
    # set up the vectors/matrices for beta estimation, availability
    for(b1 in 2:n.betas.pa){ # starts at 2 because intercept always requires a diff prior and w=1
    wpa.temp[b1] <- wpa[pos.pa[b1]]
    #wpa.temp[b1] <- 1
    mean.b.pa[b1] <- post.b.pa[b1]*(1-wpa.temp[b1])
    for(b2 in 2:n.betas.pa){ # starts at 2 because intercept always requires a diff prior and w=1
    tau.b.pa[b1,b2] <- equals(b1,b2)*((1/sd.b.pa[b1]^2)*(1-wpa.temp[b1])) + (wpa.temp[b1]*b.tau.pa)
    } # b2
    pa.beta[b1] ~ dnorm(mean.b.pa[b1],tau.b.pa[b1,b1])   # all beta coefficients
    } # b1
    
    # set up the vectors/matrices for beta estimation, perceptility
    for(b1 in 2:n.betas.pp){ # starts at 2 because intercept always requires a diff prior and w=1
    wpp.temp[b1] <- wpp[pos.pp[b1]]
    #wpp.temp[b1] <- 1
    mean.b.pp[b1] <- post.b.pp[b1]*(1-wpp.temp[b1])
    for(b2 in 2:n.betas.pp){ # starts at 2 because intercept always requires a diff prior and w=1
    tau.b.pp[b1,b2] <- equals(b1,b2)*((1/sd.b.pp[b1]^2)*(1-wpp.temp[b1])) + (wpp.temp[b1]*b.tau.pp)
    } # b2
    pp.beta[b1] ~ dnorm(mean.b.pp[b1],tau.b.pp[b1,b1])   # all beta coefficients
    } # b1
    
    # vector of stand-specific predictors
    lam.mu <- mm[,] %*% (lam.beta[]*wtemp[])
    
    stand.tau <- 1/ (stand.sig*stand.sig)
    stand.sig ~ dunif(0,10)
    
    yr.tau <- 1/ (yr.sig*yr.sig)
    yr.sig ~ dunif(0,20)
    obs.tau <- 1/ (obs.sig*obs.sig)
    obs.sig ~ dunif(0,20)
    b.tau <- 0.01
    b.tau.pa <- 0.01
    b.tau.pp <- 0.01
    
    ##### DISTANCE AND REMOVAL #####################################
    for (l in 1:L) {
    int[l] ~ dcat(pi.pa.c[site[l], yr_rot[l], ]) # removal class frequencies
    dclass[l] ~ dcat(pi.pd.c[site[l], yr_rot[l], ]) # distance class frequencies
    } # L
    
    # Distance 
    for (i in 1:nsites){
    for(t in 1:YR){  
    for(b in 1:nD){
    g[i,t,b] <- exp(-midpt[b]*midpt[b]/(2*dist.sigma[i,t]*dist.sigma[i,t])) # half-normal distance function
    f[b] <- (2*midpt[b]*delta)/(B*B)     # radial density function for point counts, change for line transects
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
    logit(p.a[i,t]) <- wpa[1]*pa.beta[1] + wpa.temp[2]*pa.beta[2]*hr[i,t] + wpa.temp[3]*pa.beta[3]*date[i,t] + 
    wpa.temp[4]*pa.beta[4]*date2[i,t] +
    wpa.temp[5]*pa.beta[5]*date[i,t]*hr[i,t] + wpa.temp[6]*pa.beta[6]*date[i,t]*date2[i,t]*hr[i,t] 
    log(dist.sigma[i,t]) <- wpp[1]*log(pp.beta[1]) +  wpp.temp[2]*pp.beta[2]*densiom[i,t] + 
    wpp.temp[3]*pp.beta[3]*noise[i,t] + wpp.temp[4]*pp.beta[4]*ba[i] + obs.eps[obs[i,t]]
    
    ##### POINT-LEVEL ABUNDANCE ###########################     
    nobs[i,t] ~ dbin(pmarg[i,t], N[i,t])  
    log(lambda[i,t]) <- lam.beta.s[stand.id[i]] + yr.eps[t] + lam.mu[i]
    N[i,t] ~ dpois(lambda[i,t])
    
    ##### GOODNESS OF FIT #######################################
    nobs.fit[i,t] ~ dbin(pmarg[i,t], N[i,t]) # create new realization of model
    e.p[i,t] <- pmarg[i,t] * N[i,t] # original model prediction
    E.p[i,t] <- pow((nobs[i,t]- e.p[i,t]),2)/(e.p[i,t]+0.5)
    E.New.p[i,t]<- pow((nobs.fit[i,t]-e.p[i,t]),2)/(e.p[i,t]+0.5)
    }} #YR #nsites 
    fit.p <- sum(E.p[1:nsites,1:YR])
    fit.new.p <- sum(E.New.p[1:nsites,1:YR])
    bayesp<-step(fit.new.p-fit.p) # Bayesian p-value for availability model. =0.5 is good fit, near 0 or 1 is poor fit

    # Random effects
    for (s in 1:S) { lam.beta.s[s] ~ dnorm(0, stand.tau) } #S
    for (y in 1:9){ yr.eps[y] ~ dnorm(0, yr.tau)}
    for (o in 1:28){ obs.eps[o] ~ dnorm(0, obs.tau)} 
    
    ##### DERIVED QUANTITIES ####################################
    for(t in 1:YR){
    Ntot[t] <- sum(N[1:nsites,t])
    D[t] <- Ntot[t] / ((3.14*B*B*nsites)/10000)  # dens per ha
    } #YR
    } # End model
    ",file="./veg-ZIP-GVS.txt")


###########################
# (5) Treatment- 
# Poisson-
# Gibb's variable selection
###########################
cat("
    model {
    ##### Variables ##########################################
    ## indices: i=site, k=visit, t=year, spp=species
    ## pa.beta = availability/removal parameters
    ## pp.beta = perceptibility/distance scale parameters
    ## dist.sigma = distance scale parameter
    ## N = detection corrected abundance
    ## Ntot = population size of total area surveyed
    ## D = density
    ## bayesp = Bayesian p-value for model fit
    
    ##### PRIORS ###############################################
    pa.beta[1] <- logit(p.pa.beta0)
    p.pa.beta0 ~ dunif(0,1) 
    pp.beta[1] ~ dunif(0, 250)
    
    # priors for the w model inclusion terms
    # this ensures that each of the 8 model combos has equal probability: Pr(m)= 1/8
    w[6] ~ dbern(0.125)
    p.w5 <- (1-w[6])*0.143 + w[6] 
    w[5] ~ dbern(p.w5) 
    p.w4 <- (1-w[6])*0.286 + w[6] 
    w[4] ~ dbern(p.w4)
    w456 <- w[4]+w[5]+w[6]
    p.w3 <- equals(w456,0)*0.5 + (1-equals(w456,0)) 
    w[3] ~ dbern(p.w3)
    w56  <- w[5]+w[6]
    p.w2 <- equals(w56,0)*0.5 + (1-equals(w56,0)) 
    w[2] ~ dbern(p.w2)
    w[1] ~ dbern(1)
    
    # priors for wpa for availability
    wpa[6] ~ dbern(0.125)
    p.wpa5 <- (1-wpa[6])*0.143 + wpa[6] 
    wpa[5] ~ dbern(p.wpa5) 
    p.wpa4 <- (1-wpa[6])*0.286 + wpa[6] 
    wpa[4] ~ dbern(p.wpa4)
    wpa456 <- wpa[4]+wpa[5]+wpa[6]
    p.wpa3 <- equals(wpa456,0)*0.5 + (1-equals(wpa456,0)) 
    wpa[3] ~ dbern(p.wpa3)
    wpa56  <- wpa[5]+wpa[6]
    p.wpa2 <- equals(wpa56,0)*0.5 + (1-equals(wpa56,0)) 
    wpa[2] ~ dbern(p.wpa2)
    wpa[1] ~ dbern(1)
    
    # priors for wpp for perceptibility
    wpp[1] ~ dbern(1) #intercept
    for (n in 2:4){wpp[n] ~ dbern(0.5)  }
    
    # set up the vectors/matrices for beta estimation, abundance
    for(b1 in 1:n.betas){
    wtemp[b1] <- w[pos[b1]]                # this uses GVS
    # wtemp[b1] <- 1                          # this forces you to fit the full model (all terms included)
    mean.b[b1] <- post.b[b1]*(1-wtemp[b1])  # prior is either 0 or full-model posterior mean
    for(b2 in 1:n.betas){                   # set up the precision matrix (inverse variance) # allows for betas to be multivariate, if desired
    tau.b[b1,b2] <- equals(b1,b2)*((1/sd.b[b1]^2)*(1-wtemp[b1])) + (wtemp[b1]*b.tau)
    } # b2
    s.beta[b1] ~ dnorm(mean.b[b1],tau.b[b1,b1])   # all beta coefficients
    } # b1
    
    # set up the vectors/matrices for beta estimation, availability
    for(b1 in 2:n.betas.pa){ # starts at 2 because intercept always requires a diff prior and w=1
    wpa.temp[b1] <- wpa[pos.pa[b1]]
    # wpa.temp[b1] <- 1
    mean.b.pa[b1] <- post.b.pa[b1]*(1-wpa.temp[b1])
    for(b2 in 2:n.betas.pa){ # starts at 2 because intercept always requires a diff prior and w=1
    tau.b.pa[b1,b2] <- equals(b1,b2)*((1/sd.b.pa[b1]^2)*(1-wpa.temp[b1])) + (wpa.temp[b1]*b.tau.pa)
    } # b2
    pa.beta[b1] ~ dnorm(mean.b.pa[b1],tau.b.pa[b1,b1])   # all beta coefficients
    } # b1
    
    # set up the vectors/matrices for beta estimation, perceptility
    for(b1 in 2:n.betas.pp){ # starts at 2 because intercept always requires a diff prior and w=1
    wpp.temp[b1] <- wpp[pos.pp[b1]]
    #  wpp.temp[b1] <- 1
    mean.b.pp[b1] <- post.b.pp[b1]*(1-wpp.temp[b1])
    for(b2 in 2:n.betas.pp){ # starts at 2 because intercept always requires a diff prior and w=1
    tau.b.pp[b1,b2] <- equals(b1,b2)*((1/sd.b.pp[b1]^2)*(1-wpp.temp[b1])) + (wpp.temp[b1]*b.tau.pp)
    } # b2
    pp.beta[b1] ~ dnorm(mean.b.pp[b1],tau.b.pp[b1,b1])   # all beta coefficients
    } # b1
    
    # vector of stand-specific predictors
    stand.mu <- mm[,] %*% (s.beta[]*wtemp[])
    
    stand.tau <- 1/ (stand.sig*stand.sig)
    stand.sig ~ dunif(0,10)
    
    yr.tau <- 1/ (yr.sig*yr.sig)
    yr.sig ~ dunif(0,20)
    obs.tau <- 1/ (obs.sig*obs.sig)
    obs.sig ~ dunif(0,20)
    b.tau <- 0.01
    b.tau.pa <- 0.01
    b.tau.pp <- 0.01
    
    ##### DISTANCE AND REMOVAL #####################################
    for (l in 1:L) {
    int[l] ~ dcat(pi.pa.c[site[l], yr_rot[l], ]) # removal class frequencies
    dclass[l] ~ dcat(pi.pd.c[site[l], yr_rot[l], ]) # distance class frequencies
    } # L
    
    # Distance 
    for (i in 1:nsites){
    for(t in 1:YR){  
    for(b in 1:nD){
    g[i,t,b] <- exp(-midpt[b]*midpt[b]/(2*dist.sigma[i,t]*dist.sigma[i,t])) # half-normal distance function
    f[b] <- (2*midpt[b]*delta)/(B*B)     # radial density function for point counts, change for line transects
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
    logit(p.a[i,t]) <- wpa[1]*pa.beta[1] + wpa.temp[2]*pa.beta[2]*hr[i,t] + wpa.temp[3]*pa.beta[3]*date[i,t] + 
    wpa.temp[4]*pa.beta[4]*date2[i,t] +
    wpa.temp[5]*pa.beta[5]*date[i,t]*hr[i,t] + wpa.temp[6]*pa.beta[6]*date[i,t]*date2[i,t]*hr[i,t] 
    log(dist.sigma[i,t]) <- wpp[1]*log(pp.beta[1]) +  wpp.temp[2]*pp.beta[2]*densiom[i,t] + 
    wpp.temp[3]*pp.beta[3]*noise[i,t] + wpp.temp[4]*pp.beta[4]*ba[i] + obs.eps[obs[i,t]]
    
    ##### POINT-LEVEL ABUNDANCE ###########################     
    nobs[i,t] ~ dbin(pmarg[i,t], N[i,t])  
    log(lambda[i,t]) <- lam.beta.s[stand.id[i]] + yr.eps[t] 
    N[i,t] ~ dpois(lambda[i,t])
    
    ##### GOODNESS OF FIT #######################################
    nobs.fit[i,t] ~ dbin(pmarg[i,t], N[i,t]) # create new realization of model
    e.p[i,t] <- pmarg[i,t] * N[i,t] # original model prediction
    E.p[i,t] <- pow((nobs[i,t]- e.p[i,t]),2)/(e.p[i,t]+0.5)
    E.New.p[i,t]<- pow((nobs.fit[i,t]-e.p[i,t]),2)/(e.p[i,t]+0.5)
    }} #YR #nsites 
    fit.p <- sum(E.p[1:nsites,1:YR])
    fit.new.p <- sum(E.New.p[1:nsites,1:YR])
    bayesp<-step(fit.new.p-fit.p) # Bayesian p-value for availability model. =0.5 is good fit, near 0 or 1 is poor fit

    # Random effects
    for (s in 1:S) { lam.beta.s[s] ~ dnorm(stand.mu[s], stand.tau) } #S
    for (y in 1:9){ yr.eps[y] ~ dnorm(0, yr.tau)}
    for (o in 1:28){ obs.eps[o] ~ dnorm(0, obs.tau)} 
    
    ##### DERIVED QUANTITIES ####################################
    for(t in 1:YR){
    Ntot[t] <- sum(N[1:nsites,t])
    D[t] <- Ntot[t] / ((3.14*B*B*nsites)/10000)  # dens per ha
    } #YR
    } # End model
    ",file="./treatment-Poisson-GVS.txt")

###########################
# (6) Treatment-
# Zero-inflated Poisson- 
# Gibb's variable selection
###########################
cat("
    model {
    ##### PRIORS ###############################################
    pa.beta[1] <- logit(p.pa.beta0)
    p.pa.beta0 ~ dunif(0,1) 
    pp.beta[1] ~ dunif(0, 250)
    
    # priors for the w model inclusion terms
    # this ensures that each of the 8 model combos has equal probability: Pr(m)= 1/8
    w[6] ~ dbern(0.125)
    p.w5 <- (1-w[6])*0.143 + w[6] 
    w[5] ~ dbern(p.w5) 
    p.w4 <- (1-w[6])*0.286 + w[6] 
    w[4] ~ dbern(p.w4)
    w456 <- w[4]+w[5]+w[6]
    p.w3 <- equals(w456,0)*0.5 + (1-equals(w456,0)) 
    w[3] ~ dbern(p.w3)
    w56  <- w[5]+w[6]
    p.w2 <- equals(w56,0)*0.5 + (1-equals(w56,0)) 
    w[2] ~ dbern(p.w2)
    w[1] ~ dbern(1)
    
    # priors for wzi for zero-inflation param
    wzi[6] ~ dbern(0.125)
    p.wzi5 <- (1-wzi[6])*0.143 + wzi[6] 
    wzi[5] ~ dbern(p.wzi5) 
    p.wzi4 <- (1-wzi[6])*0.286 + wzi[6] 
    wzi[4] ~ dbern(p.wzi4)
    wzi456 <- wzi[4]+wzi[5]+wzi[6]
    p.wzi3 <- equals(wzi456,0)*0.5 + (1-equals(wzi456,0)) 
    wzi[3] ~ dbern(p.wzi3)
    wzi56  <- wzi[5]+wzi[6]
    p.wzi2 <- equals(wzi56,0)*0.5 + (1-equals(wzi56,0)) 
    wzi[2] ~ dbern(p.wzi2)
    wzi[1] ~ dbern(1)
    
    # priors for wpa for availability
    wpa[6] ~ dbern(0.125)
    p.wpa5 <- (1-wpa[6])*0.143 + wpa[6] 
    wpa[5] ~ dbern(p.wpa5) 
    p.wpa4 <- (1-wpa[6])*0.286 + wpa[6] 
    wpa[4] ~ dbern(p.wpa4)
    wpa456 <- wpa[4]+wpa[5]+wpa[6]
    p.wpa3 <- equals(wpa456,0)*0.5 + (1-equals(wpa456,0)) 
    wpa[3] ~ dbern(p.wpa3)
    wpa56  <- wpa[5]+wpa[6]
    p.wpa2 <- equals(wpa56,0)*0.5 + (1-equals(wpa56,0)) 
    wpa[2] ~ dbern(p.wpa2)
    wpa[1] ~ dbern(1)
    
    # priors for wpp for perceptibility
    wpp[1] ~ dbern(1) #intercept
    for (n in 2:4){ wpp[n] ~ dbern(0.5)  }
    
    # set up the vectors/matrices for beta estimation, abundance
    for(b1 in 1:n.betas){
    wtemp[b1] <- w[pos[b1]]                # this uses GVS
    # wtemp[b1] <- 1                          # this forces you to fit the full model (all terms included)
    mean.b[b1] <- post.b[b1]*(1-wtemp[b1])  # prior is either 0 or full-model posterior mean
    for(b2 in 1:n.betas){                   # set up the precision matrix (inverse variance) # allows for betas to be multivariate, if desired
    tau.b[b1,b2] <- equals(b1,b2)*((1/sd.b[b1]^2)*(1-wtemp[b1])) + (wtemp[b1]*b.tau)
    } # b2
    s.beta[b1] ~ dnorm(mean.b[b1],tau.b[b1,b1])   # all beta coefficients
    } # b1
    
    # set up the vectors/matrices for beta estimation, zero-inflation
    for(b1 in 1:n.betas){
    wzi.temp[b1] <- wzi[pos[b1]]                # this uses GVS
    # wzi.temp[b1] <- 1                          # this forces you to fit the full model (all terms included)
    mean.b.zi[b1] <- post.b.zi[b1]*(1-wzi.temp[b1])  # prior is either 0 or full-model posterior mean
    for(b2 in 1:n.betas){                   # set up the precision matrix (inverse variance) # allows for betas to be multivariate, if desired
    tau.b.zi[b1,b2] <- equals(b1,b2)*((1/sd.b.zi[b1]^2)*(1-wzi.temp[b1])) + (wzi.temp[b1]*b.tau.zi)
    } # b2
    s.beta.psi[b1] ~ dnorm(mean.b.zi[b1],tau.b.zi[b1,b1])   # all beta coefficients
    } # b1
    
    # set up the vectors/matrices for beta estimation, availability
    for(b1 in 2:n.betas.pa){ # starts at 2 because intercept always requires a diff prior and w=1
    wpa.temp[b1] <- wpa[pos.pa[b1]]
    # wpa.temp[b1] <- 1
    mean.b.pa[b1] <- post.b.pa[b1]*(1-wpa.temp[b1])
    for(b2 in 2:n.betas.pa){ # starts at 2 because intercept always requires a diff prior and w=1
    tau.b.pa[b1,b2] <- equals(b1,b2)*((1/sd.b.pa[b1]^2)*(1-wpa.temp[b1])) + (wpa.temp[b1]*b.tau.pa)
    } # b2
    pa.beta[b1] ~ dnorm(mean.b.pa[b1],tau.b.pa[b1,b1])   # all beta coefficients
    } # b1
    
    # set up the vectors/matrices for beta estimation, perceptility
    for(b1 in 2:n.betas.pp){ # starts at 2 because intercept always requires a diff prior and w=1
    wpp.temp[b1] <- wpp[pos.pp[b1]]
    # wpp.temp[b1] <- 1
    mean.b.pp[b1] <- post.b.pp[b1]*(1-wpp.temp[b1])
    for(b2 in 2:n.betas.pp){ # starts at 2 because intercept always requires a diff prior and w=1
    tau.b.pp[b1,b2] <- equals(b1,b2)*((1/sd.b.pp[b1]^2)*(1-wpp.temp[b1])) + (wpp.temp[b1]*b.tau.pp)
    } # b2
    pp.beta[b1] ~ dnorm(mean.b.pp[b1],tau.b.pp[b1,b1])   # all beta coefficients
    } # b1
    
    # vector of stand-specific predictors
    stand.mu <- mm[,] %*% (s.beta[]*wtemp[]) # regressions on stand-level abundance
    stand.mu.psi <- mm[,] %*% (s.beta.psi[]*wzi.temp[]) # regressions on stand-level habitat suitability
    
    stand.tau.psi <- 1/ (stand.sig.psi*stand.sig.psi)
    stand.sig.psi ~ dunif(0,10)
    stand.tau <- 1/ (stand.sig*stand.sig)
    stand.sig ~ dunif(0,10)
    yr.tau <- 1/ (yr.sig*yr.sig)
    yr.sig ~ dunif(0,20)
    obs.tau <- 1/ (obs.sig*obs.sig)
    obs.sig ~ dunif(0,20)
    b.tau <- 0.01
    b.tau.zi <- 0.01
    b.tau.pa <- 0.01
    b.tau.pp <- 0.01
    
    ##### DISTANCE AND REMOVAL #####################################
    for (l in 1:L) {
    int[l] ~ dcat(pi.pa.c[site[l], yr_rot[l], ]) # removal class frequencies
    dclass[l] ~ dcat(pi.pd.c[site[l], yr_rot[l], ]) # distance class frequencies
    } # L
    
    # Distance 
    for (i in 1:nsites){
    for(t in 1:YR){  
    for(b in 1:nD){
    g[i,t,b] <- exp(-midpt[b]*midpt[b]/(2*dist.sigma[i,t]*dist.sigma[i,t])) # half-normal distance function
    f[b] <- (2*midpt[b]*delta)/(B*B)     # radial density function for point counts, change for line transects
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
    logit(p.a[i,t]) <- wpa[1]*pa.beta[1] + wpa.temp[2]*pa.beta[2]*hr[i,t] + wpa.temp[3]*pa.beta[3]*date[i,t] + 
    wpa.temp[4]*pa.beta[4]*date2[i,t] +
    wpa.temp[5]*pa.beta[5]*date[i,t]*hr[i,t] + wpa.temp[6]*pa.beta[6]*date[i,t]*date2[i,t]*hr[i,t] 
    log(dist.sigma[i,t]) <- wpp[1]*log(pp.beta[1]) +  wpp.temp[2]*pp.beta[2]*densiom[i,t] + wpp.temp[3]*pp.beta[3]*noise[i,t] + wpp.temp[4]*pp.beta[4]*ba[i] + obs.eps[obs[i,t]]
    
    ##### POINT-LEVEL ABUNDANCE ###########################     
    nobs[i,t] ~ dbin(pmarg[i,t], N[i,t])  
    N[i,t] ~ dpois(lambda.eff[i,t])
    lambda.eff[i,t] <- lambda[i,t] * w.lam[i,t]
    w.lam[i,t] ~ dbern(psi[i,t])
    logit(psi[i,t]) <- psi.beta.s[stand.id[i]] 
    log(lambda[i,t]) <-  lam.beta.s[stand.id[i]] + yr.eps[t]
    
    ##### GOODNESS OF FIT #######################################
    nobs.fit[i,t] ~ dbin(pmarg[i,t], N[i,t]) # create new realization of model
    e.p[i,t] <- pmarg[i,t] * N[i,t] # original model prediction
    E.p[i,t] <- pow((nobs[i,t]- e.p[i,t]),2)/(e.p[i,t]+0.5)
    E.New.p[i,t]<- pow((nobs.fit[i,t]-e.p[i,t]),2)/(e.p[i,t]+0.5)
    }} #YR #nsites 
    fit.p <- sum(E.p[1:nsites,1:YR])
    fit.new.p <- sum(E.New.p[1:nsites,1:YR])
    bayesp<-step(fit.new.p-fit.p) # Bayesian p-value for availability model. =0.5 is good fit, near 0 or 1 is poor fit

    # Random effects
    for (s in 1:S){ psi.beta.s[s] ~ dnorm(stand.mu.psi[s], stand.tau.psi)
    lam.beta.s[s] ~ dnorm(stand.mu[s], stand.tau) } #S
    for (y in 1:9){ yr.eps[y] ~ dnorm(0, yr.tau) }
    for (o in 1:28){ obs.eps[o] ~ dnorm(0, obs.tau)} # o
    
    ##### DERIVED QUANTITIES ####################################
    for(t in 1:YR){
    Ntot[t] <- sum(N[1:nsites,t])
    D[t] <- Ntot[t] / ((3.14*B*B*nsites)/10000)  # dens per ha
    } #YR
    } # End model
    " ,file="./treatment-ZIP-GVS.txt")

