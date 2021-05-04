setwd("C:\\Dropbox\\R\\Chapter2")
load("./Data/DATA.RData")
# create data frame of stand covariates
dd <- data.frame(datalfoc$CovsLam)
colnames(dd) <- c("ba2", "ba", "sf2", "sf", "dbh", "md", "scov", "scomp", "lcr")
# model matrix of stand effects (contr.sum is critical here)
mm <- model.matrix(~ 1 + ba + sf +  dbh + md + scov + scomp + lcr +
                     ba:sf + I(ba^2) + I(sf^2) + I(ba^2):sf + ba:I(sf^2),
                   data=dd)
rownames(unst)[c(7, 11, 13, 21)] <- c("md","scov", "scomp", "sf") 

# CREATE A PLOT FUNCTION THAT WILL
# AUTOMATICALLY PLOT 2 COVS FOR EACH SP
library (jagsUI)
library (scales)
plotfun<- function(sp=sp, spp=spp, cov=cov, ylims=ylims, titleline=-4,
                   ft.size=ft.size, lnwd=lnwd, ylbs=ylbs, xlbs1=xlbs1, xlbs2=xlbs2, dist="P",
                   pcol="gray60", pcol2="gray60", drive="E", ax.tcl=-0.25, ax.lwd=0.5,
                   adj.ax1=-1.5, adj.ax2=0.5, adj.ax3=1.2, adj.lab1=0.2, adj.lab3=0.2) 
{
  if (dist=="P"){
    flnm<- paste("E:\\chapter2\\outputs\\veg_GVS\\veg_global_GVS_", sp, ".Rdata", sep="")}
  else{ flnm<- paste(drive, "E:\\chapter2\\outputs\\veg_zi_GVS\\veg_global_", sp, ".Rdata", sep="") }
  #load 1st cov posteriors
  load(flnm)
  j <- which (colnames(mm)==cov[1])
  w1 <- out$sims.list$w[,j]
  beta0 <- out$sims.list$lam.beta[w1==1, 1]
  mb0 <- mean(beta0)
  beta1 <- out$sims.list$lam.beta[w1==1,j]
  mb1 <-  mean(beta1)
  # plot 1st veg cov
  ats <- (xlbs1-unst$mns[rownames(unst)==cov[1]])/unst$sds[rownames(unst)==cov[1]]
  sq1<- seq(min(mm[,cov[1]]), max(mm[,cov[1]]), length.out=100)
  pred <- array(NA, dim=c(length(sq1), length(beta1) ))
  pred.cis1 <- array(NA, dim=c(length(sq1), 3))
  for (i in 1:length(sq1)){ 
    pred[i,] <- beta0 + beta1*sq1[i] 
    pred.cis1[i ,2:3]<- quantile (pred[i,], probs = c(0.025, 0.975) )
    pred.cis1[i ,1]<- mean (pred[i,])
  }
  pred.cis1 <- exp(pred.cis1)
  plot(sq1, pred.cis1[,1], 
       xlim=c(ats[1],ats[3]), ylim=ylims, 
       lwd=lnwd, ylab="", xlab="", xaxt="n", yaxt="n", lty=1, type="n", bty="n")
  polygon(x=c(sq1, rev(sq1)), y=c(pred.cis1[,2], rev(pred.cis1[,3])), col=alpha(pcol, 0.5), border=NA) 
  lines(sq1, pred.cis1[,1], lwd=lnwd, col="black")
  #title(main=spp, line=titleline, cex.main=ft.size, font=1)
  axis(1, at=ats, labels=c(xlbs1[1], NA, xlbs1[3]), las=1, cex.axis=ft.size, cex.lab=ft.size, tcl=ax.tcl, lwd=ax.lwd, padj=adj.ax1) # c(bottom, left, top, right)
  axis(2, at=ylbs, labels=ylbs, las=1, hadj=adj.ax2, cex.axis=ft.size, cex.lab=ft.size, tcl=ax.tcl, lwd=ax.lwd)
  mtext(toupper(cov[1]), side=1, line=adj.lab1, cex=ft.size*0.66)
  box(lwd=1)
  # plot 2nd veg cov
  if(is.na(cov[2])){ 
    title(main=spp, line=titleline, cex.main=ft.size, font.main=1)
  }
  else{
    j<- which (colnames(mm)==cov[2])
    w1 <- out$sims.list$w[,j]
    beta0 <- out$sims.list$lam.beta[w1==1,1]
    mb0 <- mean(beta0)
    beta1 <- out$sims.list$lam.beta[w1==1,j]
    mb1 <-  mean(beta1)
    # plot 2nd veg cov
    ats <- (xlbs2-unst$mns[rownames(unst)==cov[2]])/unst$sds[rownames(unst)==cov[2]]
    sq<- seq(min(mm[,cov[2]]), max(mm[,cov[2]]), length.out=20)
    pred <- array(NA, dim=c(length(sq), length(beta1) ))
    pred.cis <- array(NA, dim=c(length(sq), 3))
    for (i in 1:length(sq)){ 
      pred[i,] <- beta0 + beta1*sq[i] 
      pred.cis[i ,2:3]<- quantile (pred[i,], probs = c(0.025, 0.975) )
      pred.cis[i ,1]<- mean (pred[i,])
    }
    pred.cis <- exp(pred.cis)
    par(new=T)  
    plot(sq, pred.cis[,1], 
         xlim=c(ats[1],ats[3]), ylim=ylims, 
         lwd=lnwd, ylab="", xlab="", xaxt="n", yaxt="n", lty=1, type="n", bty="n")
    polygon(x=c(sq, rev(sq)), y=c(pred.cis[,2], rev(pred.cis[,3])), col=alpha(pcol2, 0.5), border=NA) 
    lines(sq, pred.cis[,1], lwd=lnwd, col="black", lty=2)
    axis(3, at=ats, labels=c(xlbs2[1], NA, xlbs2[3]), las=1, cex.axis=ft.size, cex.lab=ft.size, tcl=ax.tcl, lwd=ax.lwd, padj=adj.ax3)
    mtext(toupper(cov[2]), side=3, line=adj.lab3, cex=ft.size*0.66)
    box(lwd=1)
    title(main=spp, line=titleline, cex.main=ft.size, font.main=1)
  }
}


# 2nd plot function to deal with quadratic basal area
plotfun.ba2<- function(sp=sp, spp=spp, cov=cov, ylims=ylims, titleline=-4,
                       ft.size=ft.size, lnwd=lnwd, ylbs=ylbs, xlbs1=xlbs1, xlbs2=xlbs2, dist="P",
                       pcol="gray60", pcol2="gray60", drive="E", ax.tcl=-0.25, ax.lwd=0.5,
                       adj.ax1=-1.5, adj.ax2=0.5, adj.ax3=1.2, adj.lab1=0.2, adj.lab3=0.2) 
{
  if (dist=="P"){
    flnm<- paste(drive, ":\\chapter2_results\\output3\\rfs_p\\", sp, "_w.Rdata", sep="")}
  else{ flnm<- paste(drive, ":\\chapter2_results\\output3\\rfs_zi\\", sp, ".Rdata", sep="") }
  #load 1st cov posteriors
  load(flnm)
  CovsLam<- as.data.frame(datalfoc$CovsLam)
  colnames(CovsLam) <- c("BA2","BA", "SPFIR2", "SPFIR", "DBH", "MID", "SCOV", "SCOMP", "LCR")
  j<- which (colnames(CovsLam)==cov[1])
  w1 <- out$sims.list$w[,j]
  w2 <- out$sims.list$w[,j-1]
  beta0 <- out$sims.list$lam.beta0[(w1*w2)==1]
  mb0 <- mean(beta0)
  beta1 <- out$sims.list$lam.beta[(w1*w2)==1,j]
  mb1 <-  mean(beta1)
  beta2 <- out$sims.list$lam.beta[(w1*w2)==1,j-1]
  # plot 1st veg cov
  ats <- (xlbs1-unst$mns[rownames(unst)==cov[1]])/unst$sds[rownames(unst)==cov[1]]
  sq1<- seq(min(CovsLam[,cov[1]])-1, max(CovsLam[,cov[1]])+1, length.out=100)
  ats2 <- (xlbs1^2-mean(veg$ba^2))/sd(veg$ba^2)
  # unstadardize the plot sequence to the natural basal area scale
  sq2 <- sq1*unst$sds[rownames(unst)==cov[1]] + unst$mns[rownames(unst)==cov[1]]
  # restandardize plot seq to the BA2 scale
  sq2 <- (sq2^2-mean(sq2^2))/sd(sq2^2)
  
  pred <- array(NA, dim=c(length(sq1), length(beta1) ))
  pred.cis1 <- array(NA, dim=c(length(sq1), 3))
  for (i in 1:length(sq1)){ pred[i,] <- beta0 + beta1*sq1[i] + beta2*sq2[i]
  pred.cis1[i ,2:3]<- quantile (pred[i,], probs = c(0.025, 0.975) )
  pred.cis1[i ,1]<- mean (pred[i,])
  }
  pred.cis1 <- exp(pred.cis1)
  plot(sq1, pred.cis1[,1], 
       xlim=c(ats[1],ats[3]), ylim=ylims, 
       lwd=lnwd, ylab="", xlab="", xaxt="n", yaxt="n", lty=1, type="n", bty="n")
  polygon(x=c(sq1, rev(sq1)), y=c(pred.cis1[,2], rev(pred.cis1[,3])), col=alpha(pcol, 0.5), border=NA) 
  lines(sq1, pred.cis1[,1], lwd=lnwd, col="black")
  #title(main=spp, line=titleline, cex.main=ft.size, font=1)
  axis(1, at=ats, labels=c(xlbs1[1], NA, xlbs1[3]), las=1, cex.axis=ft.size, cex.lab=ft.size, tcl=ax.tcl, lwd=ax.lwd, padj=adj.ax1) # c(bottom, left, top, right)
  axis(2, at=ylbs, labels=ylbs, las=1, hadj=adj.ax2, cex.axis=ft.size, cex.lab=ft.size, tcl=ax.tcl, lwd=ax.lwd)
  mtext(cov[1], side=1, line=adj.lab1, cex=ft.size*0.66)
  box(lwd=1)
  # plot 2nd veg cov
  if(is.na(cov[2])){
    title(main=spp, line=titleline, cex.main=ft.size, font.main=1)
  }
  else{
    j<- which (colnames(CovsLam)==cov[2])
    w1 <- out$sims.list$w[,j]
    beta0 <- out$sims.list$lam.beta0[w1==1]
    mb0 <- mean(beta0)
    beta1 <- out$sims.list$lam.beta[w1==1,j]
    mb1 <-  mean(beta1)
    # plot 2nd veg cov
    ats <- (xlbs2-unst$mns[rownames(unst)==cov[2]])/unst$sds[rownames(unst)==cov[2]]
    sq<- seq(min(CovsLam[,cov[2]])-1, max(CovsLam[,cov[2]])+1, length.out=20)
    pred <- array(NA, dim=c(length(sq), length(beta1) ))
    pred.cis <- array(NA, dim=c(length(sq), 3))
    for (i in 1:length(sq)){ 
      pred[i,] <- beta0 + beta1*sq[i] 
      pred.cis[i ,2:3]<- quantile (pred[i,], probs = c(0.025, 0.975) )
      pred.cis[i ,1]<- mean (pred[i,])
    }
    pred.cis <- exp(pred.cis)
    par(new=T)  
    plot(sq, pred.cis[,1], 
         xlim=c(ats[1],ats[3]), ylim=ylims, 
         lwd=lnwd, ylab="", xlab="", xaxt="n", yaxt="n", lty=1, type="n", bty="n")
    polygon(x=c(sq, rev(sq)), y=c(pred.cis[,2], rev(pred.cis[,3])), col=alpha(pcol2, 0.5), border=NA) 
    lines(sq, pred.cis[,1], lwd=lnwd, col="black", lty=2)
    axis(3, at=ats, labels=c(xlbs2[1], NA, xlbs2[3]), las=1, cex.axis=ft.size, cex.lab=ft.size, tcl=ax.tcl, lwd=ax.lwd, padj=adj.ax3)
    mtext(cov[2], side=3, line=adj.lab3, cex=ft.size*0.66)
    box(lwd=1)
    title(main=spp, line=titleline, cex.main=ft.size, font.main=1)
  }
}

# 3rd plot function to deal with quadratic SPFIR
# 2nd plot function to deal with quadratic basal area
plotfun.sf2<- function(sp=sp, spp=spp, cov=cov, ylims=ylims, titleline=-4,
                   ft.size=ft.size, lnwd=lnwd, ylbs=ylbs, xlbs1=xlbs1, xlbs2=xlbs2, dist="P",
                   pcol="gray60", pcol2="gray60", drive="E", ax.tcl=-0.25, ax.lwd=0.5,
                   adj.ax1=-1.5, adj.ax2=0.5, adj.ax3=1.2, adj.lab1=0.2, adj.lab3=0.2) 
{
  if (dist=="P"){
    flnm<- paste("E:\\chapter2\\outputs\\veg_GVS\\veg_global_GVS_", sp, ".Rdata", sep="")}
  else{ flnm<- paste(drive, "E:\\chapter2\\outputs\\veg_zi_GVS\\veg_global_", sp, ".Rdata", sep="") }
  #load 1st cov posteriors
  load(flnm)
  #j <- which (colnames(mm)==cov[1])
  w2 <- out$sims.list$w[,10]
  beta0 <- out$sims.list$lam.beta[w2==1, 1]
  mb0 <- mean(beta0)
  beta1 <- out$sims.list$lam.beta[w2==1,3]
  mb1 <-  mean(beta1)
  beta2 <- out$sims.list$lam.beta[w2==1,10]
  mb2 <-  mean(beta2)
  
  # plot 1st veg cov
  ats <- (xlbs1-unst$mns[rownames(unst)=="sf"])/unst$sds[rownames(unst)=="sf"]
  sq1<- seq(min(mm[,"sf"]), max(mm[,"sf"]), length.out=100)
  pred <- array(NA, dim=c(length(sq1), length(beta1) ))
  pred.cis1 <- array(NA, dim=c(length(sq1), 3))
  for (i in 1:length(sq1)){ 
    pred[i,] <- beta0 + beta1*sq1[i] + beta2*sq1[i]^2  
    pred.cis1[i ,2:3]<- quantile (pred[i,], probs = c(0.025, 0.975) )
    pred.cis1[i ,1]<- mean (pred[i,])
  }
  pred.cis1 <- exp(pred.cis1)
  plot(sq1, pred.cis1[,1], 
       xlim=c(ats[1],ats[3]), ylim=ylims, 
       lwd=lnwd, ylab="", xlab="", xaxt="n", yaxt="n", lty=1, type="n", bty="n")
  polygon(x=c(sq1, rev(sq1)), y=c(pred.cis1[,2], rev(pred.cis1[,3])), col=alpha(pcol, 0.5), border=NA) 
  lines(sq1, pred.cis1[,1], lwd=lnwd, col="black")
  #title(main=spp, line=titleline, cex.main=ft.size, font=1)
  axis(1, at=ats, labels=c(xlbs1[1], NA, xlbs1[3]), las=1, cex.axis=ft.size, cex.lab=ft.size, tcl=ax.tcl, lwd=ax.lwd, padj=adj.ax1) # c(bottom, left, top, right)
  axis(2, at=ylbs, labels=ylbs, las=1, hadj=adj.ax2, cex.axis=ft.size, cex.lab=ft.size, tcl=ax.tcl, lwd=ax.lwd)
  mtext(toupper("sf"), side=1, line=adj.lab1, cex=ft.size*0.66)
  box(lwd=1)
  # plot 2nd veg cov
  if(is.na(cov[2])){ 
    title(main=spp, line=titleline, cex.main=ft.size, font.main=1)
  }
  else{
    j<- which (colnames(mm)==cov[2])
    w1 <- out$sims.list$w[,j]
    beta0 <- out$sims.list$lam.beta[w1==1,1]
    mb0 <- mean(beta0)
    beta1 <- out$sims.list$lam.beta[w1==1,j]
    mb1 <-  mean(beta1)
    # plot 2nd veg cov
    ats <- (xlbs2-unst$mns[rownames(unst)==cov[2]])/unst$sds[rownames(unst)==cov[2]]
    sq<- seq(min(mm[,cov[2]]), max(mm[,cov[2]]), length.out=20)
    pred <- array(NA, dim=c(length(sq), length(beta1) ))
    pred.cis <- array(NA, dim=c(length(sq), 3))
    for (i in 1:length(sq)){ 
      pred[i,] <- beta0 + beta1*sq[i] 
      pred.cis[i ,2:3]<- quantile (pred[i,], probs = c(0.025, 0.975) )
      pred.cis[i ,1]<- mean (pred[i,])
    }
    pred.cis <- exp(pred.cis)
    par(new=T)  
    plot(sq, pred.cis[,1], 
         xlim=c(ats[1],ats[3]), ylim=ylims, 
         lwd=lnwd, ylab="", xlab="", xaxt="n", yaxt="n", lty=1, type="n", bty="n")
    polygon(x=c(sq, rev(sq)), y=c(pred.cis[,2], rev(pred.cis[,3])), col=alpha(pcol2, 0.5), border=NA) 
    lines(sq, pred.cis[,1], lwd=lnwd, col="black", lty=2)
    axis(3, at=ats, labels=c(xlbs2[1], NA, xlbs2[3]), las=1, cex.axis=ft.size, cex.lab=ft.size, tcl=ax.tcl, lwd=ax.lwd, padj=adj.ax3)
    mtext(cov[2], side=3, line=adj.lab3, cex=ft.size*0.66)
    box(lwd=1)
    title(main=spp, line=titleline, cex.main=ft.size, font.main=1)
  }
}




