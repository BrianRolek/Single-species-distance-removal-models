## ---- Prior probs veg --------
library(MuMIn)
vars <- c("ba", "sf","dbh", "md", "scov", "scomp", "lcr")
tf <- c(TRUE, FALSE)
eg <- expand.grid(tf, tf, tf, tf, tf, tf, tf)
colnames(eg) <- vars
dat <- matrix(rnorm(100*length(vars)),ncol=length(vars))
dat <- as.data.frame(dat)
names(dat) <- vars
dat$y <- rnorm(nrow(dat))

global <- lm(y ~ ba + I(ba^2) + sf + I(sf^2) + 
               dbh + md + scov + scomp + lcr +        
                ba:sf + I(ba^2):sf + ba:I(sf^2), 
                     data = dat, na.action = "na.fail")

combos <- dredge(global, eval=F, 
                 subset = 
                   dc(ba,I(ba^2),I(ba^2):sf) &&
                   dc(sf,I(sf^2),ba:I(sf^2)) &&
                   dc(ba:sf, ba:I(sf^2)) &&
                   dc(ba:sf, I(ba^2):sf)
                 )

n.mods <- length(combos)
dat2 <- dat[1,]*0 + 1
terms <- colnames(model.matrix(formula(combos[[n.mods]]),dat2))
n.terms <- length(terms)

combo.mat <- matrix(0,nrow=n.mods,ncol=n.terms)
colnames(combo.mat) <- terms

for (i in 1:n.mods){
  combo.mat[i,match(colnames(model.matrix(formula(combos[[i]]),dat2)),
                   terms)] <- 1
}
# total model probs
temp <- apply(combo.mat,2,mean) 

vegprobs <- c(w1=temp[1], w2=temp[2], w3=temp[9], w4=temp[4], w5=temp[6], 
              w6=temp[8], w7=temp[7], w8=temp[5],
              w9.ba2=mean(combo.mat[which(combo.mat[,"I(ba^2):sf"]==0),"I(ba^2)"]), 
              w10.sf2=mean(combo.mat[which(combo.mat[,"ba:I(sf^2)"]==0),"I(sf^2)"]), 
              w11.ba_sf=mean(combo.mat[which(combo.mat[,"ba:I(sf^2)"]==0 & combo.mat[,"I(ba^2):sf"]==0),"ba:sf"]), 
              w12.ba2_sf=mean(combo.mat[,"I(ba^2):sf"]) ,
              w13.ba_sf2=mean(combo.mat[,"ba:I(sf^2)"]) )
vegprobs

## ---- Prior probs treat --------
trtprobs <- matrix(c(1,0,0,0,0,0,
                    1,1,0,0,0,0,
                    1,0,1,0,0,0,
                    1,1,1,0,0,0,
                    1,0,1,1,0,0,
                    1,1,1,1,0,0,
                    1,1,1,0,1,0,
                    1,1,1,1,1,1,
                    NA, NA, NA, NA, NA, NA),
                    nrow=9,ncol=6, byrow=T, 
                    dimnames=list(c("null", "trt", "tsh", "trt+tsh", "tsh+tsh*tsh", 
                              "trt+tsh+tsh*tsh", "trt+tsh+trt*tsh", "trt+tsh+tsh*tsh+trt*tsh*tsh", "modprobs"),
                              c("w1.int","w2.trt", "w3.tsh", "w4.tsh*tsh", "w5.trt*tsh", "w6.trt*tsh*tsh")))
trtprobs[9,1] <- mean(trtprobs[,1], na.rm=T)
trtprobs[9,2] <- mean(trtprobs[1:6,2], na.rm=T)
trtprobs[9,3] <- mean(trtprobs[1:4,3], na.rm=T)
trtprobs[9,4] <- mean(trtprobs[1:7,4], na.rm=T)
trtprobs[9,5] <- mean(trtprobs[1:7,5], na.rm=T)
trtprobs[9,6] <- mean(trtprobs[1:8,6], na.rm=T)

round(trtprobs[9,], 3)

