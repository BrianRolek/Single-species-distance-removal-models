library (jagsUI)
load ("./DATA.Rdata")
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

# create data frame of stand covariates
dd <- data.frame(datalfoc$CovsLam)
colnames(dd) <- c("ba2", "ba", "sf2", "sf", "dbh", "md", "scov", "scomp", "lcr")
# model matrix of stand effects (contr.sum is critical here)
mmzi <- mm <- model.matrix(~ 1 + ba + sf +  dbh + md + scov + scomp + lcr +
                     ba:sf + I(ba^2) + I(sf^2) + I(ba^2):sf + ba:I(sf^2),
                     data=dd)

# position of the beta coefficients associated with each bernoilli indicator var (i.e., treat has 7 terms w/ intercept)
pos.zi <- pos <- as.numeric(attr(mm,"assign")+1)
n.betas.zi <- n.betas <- length(pos)
pos.pa <- c(1:6)
n.betas.pa <- length(pos.pa)
pos.pp <- c(1:4)
n.betas.pp <- length(pos.pp)

# Define model in BUGS
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
    psi.beta[1] <- logit(p.psi.beta0)
    p.psi.beta0 ~ dunif(0,1)
    
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
    
    # priors for the w model inclusion terms
    wzi[13] ~ dbern(0.1667)
    wzi[12] ~ dbern(0.1667)
    wzi12_13 <- wzi[12]+wzi[13]
    p.wzi11 <- equals(wzi12_13,0)*.3077 + (1-equals(wzi12_13,0))
    wzi[11] ~ dbern(p.wzi11)
    p.wzi10 <- equals(wzi[13],0)*0.4 + wzi[13]
    wzi[10] ~ dbern(p.wzi10)
    p.wzi9 <- equals(wzi[12],0)*0.4 + wzi[12]
    wzi[9] ~ dbern(p.wzi9)
    wzi[8] ~ dbern(0.5)
    wzi[7] ~ dbern(0.5)
    wzi[6] ~ dbern(0.5)
    wzi[5] ~ dbern(0.5)
    wzi[4] ~ dbern(0.5)
    wzi10_13 <- wzi[10] + wzi[11] + wzi[12] + wzi[13]
    p.wzi3 <- equals(wzi10_13,0)*0.5 + (1-equals(wzi10_13,0))
    wzi[3] ~ dbern(p.wzi3)
    wzi9.11_13 <- wzi[9] + wzi[11] + wzi[12] + wzi[13]
    p.wzi2 <- equals(wzi9.11_13,0)*0.5 + (1-equals(wzi9.11_13,0))
    wzi[2] ~ dbern(p.wzi2)
    wzi[1] ~ dbern(1)
    wzi.temp[1] <- 1

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
    #wtemp[b1] <- w[pos[b1]]                # this uses GVS
    wtemp[b1] <- 1                          # this forces you to fit the full model (all terms included)
    mean.b[b1] <- post.b[b1]*(1-wtemp[b1])  # prior is either 0 or full-model posterior mean
    for(b2 in 1:n.betas){                   # set up the precision matrix (inverse variance) # allows for betas to be multivariate, if desired
    tau.b[b1,b2] <- equals(b1,b2)*((1/sd.b[b1]^2)*(1-wtemp[b1])) + (wtemp[b1]*b.tau)
    } # b2
    lam.beta[b1] ~ dnorm(mean.b[b1],tau.b[b1,b1])   # all beta coefficients
    } # b1
    
    # set up the vectors/matrices for beta estimation, abundance
    for(b1 in 2:n.betas.zi){  # starts at 2 because intercept always requires a diff prior and w=1
    #wzi.temp[b1] <- wzi[pos.zi[b1]]                # this uses GVS
    wzi.temp[b1] <- 1                          # this forces you to fit the full model (all terms included)
    mean.b.zi[b1] <- post.b.zi[b1]*(1-wzi.temp[b1])  # prior is either 0 or full-model posterior mean
    for(b2 in 2:n.betas.zi){                   # set up the precision matrix (inverse variance) # allows for betas to be multivariate, if desired
    tau.b.zi[b1,b2] <- equals(b1,b2)*((1/sd.b.zi[b1]^2)*(1-wzi.temp[b1])) + (wzi.temp[b1]*b.tau.zi)
    } # b2
    psi.beta[b1] ~ dnorm(mean.b.zi[b1],tau.b.zi[b1,b1])   # all beta coefficients
    } # b1

    # set up the vectors/matrices for beta estimation, availability
    for(b1 in 2:n.betas.pa){ # starts at 2 because intercept always requires a diff prior and w=1
    #wpa.temp[b1] <- wpa[pos.pa[b1]]
    wpa.temp[b1] <- 1
    mean.b.pa[b1] <- post.b.pa[b1]*(1-wpa.temp[b1])
    for(b2 in 2:n.betas.pa){ # starts at 2 because intercept always requires a diff prior and w=1
    tau.b.pa[b1,b2] <- equals(b1,b2)*((1/sd.b.pa[b1]^2)*(1-wpa.temp[b1])) + (wpa.temp[b1]*b.tau.pa)
    } # b2
    pa.beta[b1] ~ dnorm(mean.b.pa[b1],tau.b.pa[b1,b1])   # all beta coefficients
    } # b1
    
    # set up the vectors/matrices for beta estimation, perceptility
    for(b1 in 2:n.betas.pp){ # starts at 2 because intercept always requires a diff prior and w=1
    #wpp.temp[b1] <- wpp[pos.pp[b1]]
    wpp.temp[b1] <- 1
    mean.b.pp[b1] <- post.b.pp[b1]*(1-wpp.temp[b1])
    for(b2 in 2:n.betas.pp){ # starts at 2 because intercept always requires a diff prior and w=1
    tau.b.pp[b1,b2] <- equals(b1,b2)*((1/sd.b.pp[b1]^2)*(1-wpp.temp[b1])) + (wpp.temp[b1]*b.tau.pp)
    } # b2
    pp.beta[b1] ~ dnorm(mean.b.pp[b1],tau.b.pp[b1,b1])   # all beta coefficients
    } # b1
    
    # vector of site-specific predictors
    lam.mu <- mm[,] %*% (lam.beta[]*wtemp[])
    psi.mu <- mmzi[,] %*% (psi.beta[]*wzi.temp[])
    
    stand.tau <- 1/ (stand.sig*stand.sig)
    stand.sig ~ dunif(0,10)
    psi.stand.tau <- 1/ (psi.stand.sig*psi.stand.sig)
    psi.stand.sig ~ dunif(0,10)
    yr.tau <- 1/ (yr.sig*yr.sig)
    yr.sig ~ dunif(0,20)
    psi.yr.tau <- 1/ (psi.yr.sig*psi.yr.sig)
    psi.yr.sig ~ dunif(0,20)
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
    logit(p.a[i,t]) <- wpa[1]*pa.beta[1] + wpa.temp[2]*pa.beta[2]*hr[i,t] + wpa.temp[3]*pa.beta[3]*date[i,t] + 
                        wpa.temp[4]*pa.beta[4]*date2[i,t] +
                         wpa.temp[5]*pa.beta[5]*date[i,t]*hr[i,t] + wpa.temp[6]*pa.beta[6]*date[i,t]*date2[i,t]*hr[i,t] 
    log(dist.sigma[i,t]) <- wpp[1]*log(pp.beta[1]) +  wpp.temp[2]*pp.beta[2]*densiom[i,t] + 
                            wpp.temp[3]*pp.beta[3]*noise[i,t] + wpp.temp[4]*pp.beta[4]*ba[i] + obs.eps[obs[i,t]]
    
    ##### POINT-LEVEL ABUNDANCE ###########################     
    nobs[i,t] ~ dbin(pmarg[i,t], N[i,t])  
    N[i,t] ~ dpois(lambda.eff[i,t])
    lambda.eff[i,t] <- lambda[i,t] * w.lam[i,t]
    w.lam[i,t] ~  dbern(psi[i,t])
    log(lambda[i,t]) <- lam.beta.s[stand.id[i]] + yr.eps[t] + lam.mu[i]
    logit(psi[i,t]) <- psi.mu[i]
    
    ##### GOODNESS OF FIT #######################################
    nobs.fit[i,t] ~ dbin(pmarg[i,t], N[i,t]) # create new realization of model
    e.p[i,t] <- pmarg[i,t] * N[i,t] # original model prediction
    E.p[i,t] <- pow((nobs[i,t]- e.p[i,t]),2)/(e.p[i,t]+0.5)
    E.New.p[i,t]<- pow((nobs.fit[i,t]-e.p[i,t]),2)/(e.p[i,t]+0.5)
    }} #YR #nsites 
    
    for (s in 1:S) {
    lam.beta.s[s] ~ dnorm(0, stand.tau)
    } #S
    # Random effects
    for (y in 1:9){ 
      yr.eps[y] ~ dnorm(0, yr.tau)
      } # y
    for (o in 1:28){ obs.eps[o] ~ dnorm(0, obs.tau)} 
    
    ##### DERIVED QUANTITIES ####################################
    for(t in 1:YR){
    Ntot[t] <- sum(N[1:nsites,t])
    D[t] <- Ntot[t] / ((3.14*B*B*nsites)/10000)  # dens per ha
    } #YR
    
    fit.p <- sum(E.p[1:nsites,1:YR])
    fit.new.p <- sum(E.New.p[1:nsites,1:YR])
    bayesp <- step(fit.new.p-fit.p) # Bayesian p-value for availability model. =0.5 is good fit, near 0 or 1 is poor fit
    } # End model
    " ,file="./model_V-zip_global.txt")

for (i in c(11,10,3,5,2)){ #use zips for spp c(11,10,3,5,2)
  try(rm("out"))
  spp <- spp.list.foc[i]
  spp.num<- which(dimnames(nobs)[[3]]==spp)
  datalfoc$nobs <- Nav <- apply(ab2[,1:2,,,spp], c(1,4),sum, na.rm=T)
  Mst <- apply(Nav, c(1), max, na.rm=T) +1
  
  inits <- function(){  list(
    N = Nav,
    p.psi.beta0= runif(1, 0.9, 1),
    p.pa.beta0= runif(1, 0.3, 0.8),
    pp.beta= c(runif(1, 20, 65), rep(NA, datalfoc$nCovsPP-2)),
    stand.sig= runif(1, 0, 2), 
    s.beta = runif(n.betas,-.5,.5)
  )  }
  
  params <- c("pa.beta", "pp.beta", "obs.eps", "obs.sig",
              "lam.beta", "stand.sig", "yr.eps", "yr.sig",  
              "psi.beta", 
              "Ntot", "D", 
              "bayesp", "w", "wzi", "wpa", "wpp"
  )
  
  datalfoc$nobs <- nobs[,,spp]
  datalfoc$dclass <- dclass[datalfoc$species==spp.num]
  datalfoc$int <- int[datalfoc$species==spp.num]
  datalfoc$site <- site[datalfoc$species==spp.num]
  datalfoc$yr_rot <- yr_rot[datalfoc$species==spp.num]
  datalfoc$L <- length(datalfoc$species[datalfoc$species==spp.num])
  
  # new objects for beta estimation
  datalfoc$mm <- mm
  datalfoc$mmzi <- mmzi
  datalfoc$pos <- pos
  datalfoc$pos.zi <- pos.zi
  datalfoc$pos.pa <- pos.pa
  datalfoc$pos.pp <- pos.pp
  datalfoc$n.betas <- n.betas
  datalfoc$n.betas.zi <- n.betas.zi
  datalfoc$n.betas.pa <- n.betas.pa
  datalfoc$n.betas.pp <- n.betas.pp

  datalfoc$post.b <- rep(0, n.betas) # out$mean$post.b
  datalfoc$sd.b <- rep(10, n.betas) 
  datalfoc$post.b.zi <- rep(0, n.betas.zi) # out$mean$post.b
  datalfoc$sd.b.zi <- rep(10, n.betas.zi) 
  datalfoc$post.b.pa <- rep(0, n.betas.pa) # out$mean$post.b
  datalfoc$sd.b.pa <- rep(10, n.betas.pa) 
  datalfoc$post.b.pp <- rep(0, n.betas.pp) # out$mean$post.b
  datalfoc$sd.b.pp <- rep(10, n.betas.pp) 
  # MCMC settings
  ni <- 200000  ;   nb <- 100000   ;   nt <- 10   ;   nc <- 3 ; na=10000
#  ni <- 100 ;   nb <- 50   ;   nt <- 1   ;   nc <- 1 ; na <- 100
  # Run JAGS
  out <- jags(datalfoc, inits=inits, 
              params, "./model_V-zip_global.txt",
             # "V14.txt", 
              n.thin=nt, n.chains=nc, 
              n.burnin=nb, n.iter=ni, n.adapt=na,
              parallel = T, modules=c("glm"))
  
  fn<- paste( "./", spp, "_V-zip_global.RData", sep="" )
  save(list= c("out", "datalfoc"), file=fn)
}

