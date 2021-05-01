############################
# Supplemental materials for:
# B. W. Rolek, D. J. Harrison, D. W. Linden,  C. S. Loftin, 
# P. B. Wood. 2021. Associations among breeding 
# conifer-associated birds, forestry treatments, 
# years-since-harvest, and vegetation characteristics in 
# regenerating stands. 
#############################
# software used
# JAGS 4.3.0 
# R version 4.0.2 (2020-06-22) -- "Taking Off Again"
## ---- treatment Poisson GVS --------
library (jagsUI) # v1.5.1
load("./global_est.Rdata") # output file from global model
load ("./DATA.RData")
# data manipulation
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
dd <- data.frame(treat=factor(datalfoc$treat),tsh=datalfoc$tsh,tsh2=datalfoc$tsh^2)
# model matrix of stand effects (contr.sum is critical here)
mm <- model.matrix(~treat*tsh+treat*tsh2,dd,contrasts=list(treat="contr.sum"))
# position of the beta coefficients associated with bernoulli indicator variable (i.e., treat has 7 terms w/ intercept)
pos <- as.numeric(attr(mm,"assign")+1)
n.betas <- length(pos)
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
      log(lambda[i,t]) <- lam.beta.s[stand.id[i]] + yr.eps[t] 
      N[i,t] ~ dpois(lambda[i,t])

##### GOODNESS OF FIT #######################################
nobs.fit[i,t] ~ dbin(pmarg[i,t], N[i,t]) # create new realization of model
e.p[i,t] <- pmarg[i,t] * N[i,t] # original model prediction
E.p[i,t] <- pow((nobs[i,t]- e.p[i,t]),2)/(e.p[i,t]+0.5)
E.New.p[i,t]<- pow((nobs.fit[i,t]-e.p[i,t]),2)/(e.p[i,t]+0.5)
    }} #YR #nsites 

for (s in 1:S) {
  lam.beta.s[s] ~ dnorm(stand.mu[s], stand.tau)
} #S
# Random effects
for (y in 1:9){ yr.eps[y] ~ dnorm(0, yr.tau)}
for (o in 1:28){ obs.eps[o] ~ dnorm(0, obs.tau)} 

##### DERIVED QUANTITIES ####################################
for(t in 1:YR){
      Ntot[t] <- sum(N[1:nsites,t])
      D[t] <- Ntot[t] / ((3.14*B*B*nsites)/10000)  # dens per ha
      } #YR

fit.p <- sum(E.p[1:nsites,1:YR])
fit.new.p <- sum(E.New.p[1:nsites,1:YR])
bayesp<-step(fit.new.p-fit.p) # Bayesian p-value for availability model. =0.5 is good fit, near 0 or 1 is poor fit
    } # End model
    "
    ,file="./T-pois-GVS.txt")
# CAUTION: These next lines run the model and
# take a VERY long time to run 
# (>1 week, each species took 4 days)
# We ran these on an HPC and specified the
# loop to run 4-5 species sequentially. 
for (i in 1:19){ #Create 5 files: 1:4, 5:8, 9:12, 13:16, 17:19
try(rm("out"))
spp <- spp.list.foc[i]
spp.num<- which(dimnames(nobs)[[3]]==spp)
datalfoc$nobs <- Nav <- apply(ab2[,1:2,,,spp], c(1,4),sum, na.rm=T)
Mst <- apply(Nav, c(1), max, na.rm=T) +1

inits <- function(){  list(
  N = Nav,
  p.pa.beta0= runif(1, 0.3, 0.8),
  pp.beta= c(runif(1, 20, 65), rep(NA, datalfoc$nCovsPP-2)),
  stand.sig= runif(1, 0, 2), 
  s.beta = runif(n.betas,-.5,.5)
  )  }

params <- c("pa.beta", "pp.beta", 
            "lam.beta", "lam.beta1", "lam.beta2", 
            "Ntot", "D", 
            "stand.sig", "s.beta", 
            "bayesp", "w", "wpa", "wpp",
            "yr.eps", "yr.sig", "obs.eps", "obs.sig"
            )

datalfoc$nobs <- nobs[,,spp]
datalfoc$dclass <- dclass[datalfoc$species==spp.num]
datalfoc$int <- int[datalfoc$species==spp.num]
datalfoc$site <- site[datalfoc$species==spp.num]
datalfoc$yr_rot <- yr_rot[datalfoc$species==spp.num]
datalfoc$L <- length(datalfoc$species[datalfoc$species==spp.num])

# new objects for beta estimation
datalfoc$mm <- mm
datalfoc$pos <- pos
datalfoc$pos.pa <- pos.pa
datalfoc$pos.pp <- pos.pp
datalfoc$n.betas <- n.betas
datalfoc$n.betas.pa <- n.betas.pa
datalfoc$n.betas.pp <- n.betas.pp
# input posterior means & sds from a global model run when using GVS
spp.num2 <- which( names(post.b)==spp )
datalfoc$post.b <- post.b[[spp.num2]] # out$mean$post.b
datalfoc$sd.b <-sd.b[[spp.num2]] 
datalfoc$post.b.pa <- post.b.pa[[spp.num2]] # out$mean$post.b
datalfoc$sd.b.pa <- sd.b.pa[[spp.num2]]
datalfoc$post.b.pp <- post.b.pp[[spp.num2]] # out$mean$post.b
datalfoc$sd.b.pp <- sd.b.pp[[spp.num2]]

# MCMC settings
ni <- 200000  ;   nb <- 100000   ;   nt <- 10   ;   nc <- 6 ; na=10000
#ni <- 100 ;   nb <- 50   ;   nt <- 1   ;   nc <- 1 ; na <- 100
# Run JAGS
out <- jags(datalfoc, inits=inits, 
            params, "./T-pois-GVS.txt",
            n.thin=nt, n.chains=nc, 
            n.burnin=nb, n.iter=ni, n.adapt=na
            , parallel = T, modules=c("glm"),
            codaOnly= "N")

fn<- paste( "./", spp, "_T-pois-GVS.RData", sep="" )
save(list= c("out", "datalfoc"), file=fn)
}

