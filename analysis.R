library (jagsUI)
load (".\\DATA.Rdata")
load(".\\global_est_veg.Rdata")
source(".\\models.R")
datalfoc$SPP <- length(spp.list.foc)
yr <- array(NA, dim=c(dim (ab)[1], 9) )
yr[,1:3] <- 1; yr[,4:6] <- 2; yr[,7:9] <- 3
datalfoc$yr <- yr
datalfoc$tsh.pred <- seq(length.out=7, -3, 3)
datalfoc$tsh2.pred <- seq(length.out=7, -3, 3)
s.year <- array(NA, dim=c(114, 9))
s.year[,1:3] <- 1; s.year[,4:6] <- 2; s.year[,7:9] <- 3
datalfoc$s.year <- s.year

datalfoc$ba <- datalfoc$CovsLam[, "ba"]
nobs <- datalfoc$nobs
dclass <- datalfoc$dclass
int <- datalfoc$int
site <- datalfoc$site
yr_rot <- datalfoc$yr_rot

# print number of detections
apply(ab2[,1:2,,,dimnames(ab2)[[5]] %in% spp.list.foc], c(5), sum, na.rm=T)

# create data frame of stand covariates
dd <- data.frame(datalfoc$CovsLam)
colnames(dd) <- c("ba2", "ba", "sf2", "sf", "dbh", "md", "scov", "scomp", "lcr")
# model design matrix of covariates
mm <- model.matrix(~ 1 + ba + sf +  dbh + md + scov + scomp + lcr +
                     ba:sf + I(ba^2) + I(sf^2) + I(ba^2):sf + ba:I(sf^2),
                   data=dd)

# position of the beta coefficients associated with each term (i.e., treat has 7 terms w/ intercept)
pos <- as.numeric(attr(mm,"assign")+1)
n.betas <- length(pos)
pos.pa <- c(1:6)
n.betas.pa <- length(pos.pa)
pos.pp <- c(1:4)
n.betas.pp <- length(pos.pp)

################################
# This will take a very long time to run
# 4 days per species on a cluster
# Or break into smaller runs
################################
for (i in 1:19){ #Skip spp c(10,3,5,2) use ZIPs
  try(rm("out"))
  spp <- spp.list.foc[i]
  spp.num<- which(dimnames(nobs)[[3]]==spp)
  # Inits and parameters to save
  # Crunch the numbers, reformat
  datalfoc$nobs <- Nav <- apply(ab2[,1:2,,,spp], c(1,4),sum, na.rm=T)
  Mst <- apply(Nav, c(1), max, na.rm=T) +1
  
  # remove comments to add covariates
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
              "stand.sig", "lam.beta", 
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
  # these should be replaced with actual posterior means & sds from a full model run!!
  spp.num2 <- which( names(post.b)==spp )
  datalfoc$post.b <- post.b[[spp.num2]] # out$mean$post.b
  datalfoc$sd.b <-sd.b[[spp.num2]]
  datalfoc$post.b.pa <- post.b.pa[[spp.num2]] # out$mean$post.b
  datalfoc$sd.b.pa <- sd.b.pa[[spp.num2]]
  datalfoc$post.b.pp <- post.b.pp[[spp.num2]] # out$mean$post.b
  datalfoc$sd.b.pp <- sd.b.pp[[spp.num2]]
  
  # uncomment when running run global model
  # datalfoc$post.b <- rep(0, n.betas) # out$mean$post.b
  # datalfoc$sd.b <- rep(10, n.betas) 
  # datalfoc$post.b.pa <- rep(0, n.betas.pa) # out$mean$post.b
  # datalfoc$sd.b.pa <- rep(10, n.betas.pa) 
  # datalfoc$post.b.pp <- rep(0, n.betas.pp) # out$mean$post.b
  # datalfoc$sd.b.pp <- rep(10, n.betas.pp) 
  
  # MCMC settings
  ni <- 200000  ;   nb <- 100000   ;   nt <- 10   ;   nc <- 6 ; na=10000
  # Run JAGS
  out <- jags(datalfoc, inits=inits, 
              params, "./veg-Poisson-GVS.txt",
              #"V14.txt", 
              n.thin=nt, n.chains=nc, 
              n.burnin=nb, n.iter=ni, n.adapt=na,
              parallel = T, modules=c("glm"),
              codaOnly= "N")
  
  fn<- paste( "./veg-Poisson-GVS_", spp, ".RData", sep="" )
  save(list= c("out", "datalfoc"), file=fn)
}


