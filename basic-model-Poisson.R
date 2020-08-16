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
    f[i,t,b] <- (2*midpt[b]*delta)/(B*B)     # radial density function for point counts, change for line transects
    pi.pd[i,t,b] <- g[i,t,b]*f[i,t,b]
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
    
    ##### DERIVED QUANTITIES ####################################
    for(t in 1:YR){
    Ntot[t] <- sum(N[1:nsites,t])
    D[t] <- Ntot[t] / ((3.14*B*B*nsites)/10000)  # dens per ha
    } #YR
    
    fit.p <- sum(E.p[1:nsites,1:YR])
    fit.new.p <- sum(E.New.p[1:nsites,1:YR])
    bayesp<-step(fit.new.p-fit.p) # Bayesian p-value for availability model. =0.5 is good fit, near 0 or 1 is poor fit
    } # End model
    ",file="./basic-model-Poisson.txt")