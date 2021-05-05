###########################
# each panel as a species
source("C:\\Users\\rolek.brian\\Documents\\Projects\\CH2\\Results\\plot-functions.R")
fname<- "C:\\Users\\rolek.brian\\Documents\\Projects\\CH2\\docs\\veg_bySpecies.pdf"
pdf(fname, height=6, width=6) # width 36
par (mfrow=c(5,4), mar=c(1.25, 2, 1.25, 0), oma=c(3.5,3.5,0.5,2)) # c(bottom, left, top, right)
ft.size <- 1
lnwd <- 3
ttl <- -2
pcol="blue"
pcol2="red"

plotfun(sp="BBWA", cov=c("scov", "sf"), spp="Bay-breasted\nWarbler", titleline=ttl, pcol=pcol, pcol2=pcol2,
        ylims=c(0,0.2), ft.size=ft.size, lnwd=lnwd, ylbs=c(0,0.1,0.2), xlbs1=c(0,0.45,0.9), xlbs2=c(0,0.5,1),
        adj.ax1=-1.5, adj.ax2=0.5, adj.ax3=1.2, adj.lab1=0.2, adj.lab3=0.2, drive=drive) 

plotfun(sp="BBWA", cov=c("lcr", "md"), spp="Bay-breasted\nWarbler", titleline=ttl, pcol=pcol, pcol2=pcol2,
        ylims=c(0,0.2), ft.size=ft.size, lnwd=lnwd, ylbs=c(0,0.1,0.2), xlbs1=c(0,0.5,1), xlbs2=c(0,0.5,1),
        adj.ax1=-1.5, adj.ax2=0.5, adj.ax3=1.2, adj.lab1=0.2, adj.lab3=0.2, drive=drive)

plotfun(sp="BLBW", cov=c("dbh", "ba"), spp="Blackburnian\nWarbler", titleline=ttl, pcol=pcol, pcol2=pcol2,
        ylims=c(0,3), ft.size=ft.size, lnwd=lnwd, ylbs=c(0,1.5,3), xlbs1=c(0,30,60), xlbs2=c(0,35,70),
        adj.ax1=-1.5, adj.ax2=0.5, adj.ax3=1.2, adj.lab1=0.2, adj.lab3=0.2, drive=drive)

plotfun(sp="BLBW", cov=c("md"), spp="Blackburnian\nWarbler", titleline=ttl, pcol=pcol, pcol2=pcol2,
        ylims=c(0,0.8), ft.size=ft.size, lnwd=lnwd, ylbs=c(0,0.4,0.8), xlbs1=c(0,0.5,1),
        adj.ax1=-1.5, adj.ax2=0.5, adj.ax3=1.2, adj.lab1=0.2, adj.lab3=0.2, drive=drive) 

plotfun(sp="BLPW", cov=c("ba", "lcr"), spp="Blackpoll\nWarbler", titleline=ttl, pcol=pcol, pcol2=pcol2,
        ylims=c(0,0.3), ft.size=ft.size, lnwd=lnwd, ylbs=c(0,0.15,0.3), xlbs1=c(0,35,70), xlbs2=c(0,0.5,1),
        adj.ax1=-1.5, adj.ax2=0.5, adj.ax3=1.2, adj.lab1=0.2, adj.lab3=0.2, drive=drive) 

plotfun(sp="CAWA", cov=c("scomp","ba"), spp="Canada\nWarbler", titleline=ttl, pcol=pcol, pcol2=pcol2,
        ylims=c(0,1.0), ft.size=ft.size, lnwd=lnwd, ylbs=c(0,0.5,1.0), xlbs1=c(0,0.45,0.9), xlbs2=c(0,35,70),
        adj.ax1=-1.5, adj.ax2=0.5, adj.ax3=1.2, adj.lab1=0.2, adj.lab3=0.2, drive=drive)  

plotfun(sp="HETH", cov=c("sf"), spp="Hermit\nThrush", titleline=ttl, pcol=pcol, pcol2=pcol2,
        ylims=c(0,0.8), ft.size=ft.size, lnwd=lnwd, ylbs=c(0,0.4,0.8), xlbs1=c(0,0.5,1),
        adj.ax1=-1.5, adj.ax2=0.5, adj.ax3=1.2, adj.lab1=0.2, adj.lab3=0.2, drive=drive) 

plotfun(sp="OSFL", cov=c("md"), spp="Olive-sided\nFlycatcher", titleline=ttl, pcol=pcol, pcol2=pcol2,
        ylims=c(0,0.03), ft.size=ft.size, lnwd=lnwd, ylbs=c(0,0.01,0.02,0.03), xlbs1=c(0,0.5,1),
        adj.ax1=-1.5, adj.ax2=0.5, adj.ax3=1.2, adj.lab1=0.2, adj.lab3=0.2, drive=drive) 

plotfun(sp="YPWA", cov=c("ba", "md"), spp="Palm\nWarbler", titleline=ttl, pcol=pcol, pcol2=pcol2,
        ylims=c(0,0.5), ft.size=ft.size, lnwd=lnwd, ylbs=c(0,0.25,0.5), xlbs1=c(0,35,70), xlbs2=c(0,0.5,1.0),
        adj.ax1=-1.5, adj.ax2=0.5, adj.ax3=1.2, adj.lab1=0.2, adj.lab3=0.2, drive=drive)

plotfun(sp="RBNU", cov=c("dbh","ba"), spp="Red-breasted\nNuthatch", titleline=ttl, pcol=pcol, pcol2=pcol2,
        ylims=c(0,1.5), ft.size=ft.size, lnwd=lnwd, ylbs=c(0,0.75,1.5), xlbs1=c(0,30,60), xlbs2=c(0,35,70),
        adj.ax1=-1.5, adj.ax2=0.5, adj.ax3=1.2, adj.lab1=0.2, adj.lab3=0.2, drive=drive)

plotfun.sf2(sp="RBNU", cov=c("sf2", NA), spp="Red-breasted\nNuthatch", titleline=ttl, pcol=pcol, pcol2=pcol2,
            ylims=c(0,0.6), ft.size=ft.size, lnwd=lnwd, ylbs=c(0,0.3,0.6), xlbs1=c(0,0.5,1), 
            adj.ax1=-1.5, adj.ax2=0.5, adj.ax3=1.2, adj.lab1=0.2, adj.lab3=0.2, drive=drive)

plotfun(sp="RCKI", cov=c("sf","ba"), spp="Ruby-crowned\nKinglet", titleline=ttl, pcol=pcol, pcol2=pcol2,
        ylims=c(0,0.5), ft.size=ft.size, lnwd=lnwd, ylbs=c(0,0.25,0.5), xlbs1=c(0,0.5,1), xlbs2=c(0,35,70),
        adj.ax1=-1.5, adj.ax2=0.5, adj.ax3=1.2, adj.lab1=0.2, adj.lab3=0.2, drive=drive)

plotfun(sp="SWTH", cov=c("sf"), spp="Swainson's\nThrush", titleline=ttl, pcol=pcol, pcol2=pcol2,
        ylims=c(0,0.8), ft.size=ft.size, lnwd=lnwd, ylbs=c(0,0.4,0.8), xlbs1=c(0,0.5,1),
        adj.ax1=-1.5, adj.ax2=0.5, adj.ax3=1.2, adj.lab1=0.2, adj.lab3=0.2, drive=drive) 

plotfun(sp="WIWR", cov=c("sf"), spp="Winter\nWren", titleline=ttl, pcol=pcol, pcol2=pcol2,
        ylims=c(0,0.6), ft.size=ft.size, lnwd=lnwd, ylbs=c(0,0.3,0.6), xlbs1=c(0,0.5,1), 
        adj.ax1=-1.5, adj.ax2=0.5, adj.ax3=1.2, adj.lab1=0.2, adj.lab3=0.2, drive=drive)

plotfun(sp="WTSP", cov=c("scov","md"), spp="White-throated\nSparrow", titleline=ttl, pcol=pcol, pcol2=pcol2,
        ylims=c(0,2), ft.size=ft.size, lnwd=lnwd, ylbs=c(0,1,2), xlbs1=c(0,0.45,0.9), xlbs2=c(0,0.5,1),
        adj.ax1=-1.5, adj.ax2=0.5, adj.ax3=1.2, adj.lab1=0.2, adj.lab3=0.2, drive=drive) 

plotfun.sf2(sp="WTSP", cov=c("sf2", NA), spp="White-throated\nSparrow", titleline=ttl, pcol=pcol, pcol2=pcol2,
            ylims=c(0,0.8), ft.size=ft.size, lnwd=lnwd, ylbs=c(0,0.4,0.8), xlbs1=c(0,0.5,1), 
            adj.ax1=-1.5, adj.ax2=0.5, adj.ax3=1.2, adj.lab1=0.2, adj.lab3=0.2, drive=drive)

plotfun(sp="MYWA", cov=c("scov","sf"), spp="Yellow-rumped\nWarbler", titleline=ttl, pcol=pcol, pcol2=pcol2,
        ylims=c(0,1), ft.size=ft.size, lnwd=lnwd, ylbs=c(0,0.5,1), xlbs1=c(0,0.45,0.9), xlbs2=c(0,0.5,1),
        adj.ax1=-1.5, adj.ax2=0.5, adj.ax3=1.2, adj.lab1=0.2, adj.lab3=0.2, drive=drive)

plotfun(sp="MYWA", cov=c("lcr"), spp="Yellow-rumped\nWarbler", titleline=ttl, pcol=pcol, pcol2=pcol2,
        ylims=c(0,1), ft.size=ft.size, lnwd=lnwd, ylbs=c(0,0.5,1.0), xlbs1=c(0,0.5,1),
        adj.ax1=-1.5, adj.ax2=0.5, adj.ax3=1.2, adj.lab1=0.2, adj.lab3=0.2, drive=drive)

plot(NA, new=T)
legend(0, 1.0, legend=c("Upper x-axis", "Lower x-axis"), 
       fill=c(pcol2, pcol), border=NA, 
       xpd=NA, lwd=3, lty=c(2,1), 
       cex=1.5, y.intersp=1.5, bty="n")
mtext("Abundance", side=2, line=1, cex=1, outer=T)
mtext("Vegetation covariates", side=1, line=1, cex=1, outer=T)
dev.off()

##############
# Coefficient plots
##############
library (viridis)
cdat <- read.csv("C:\\Users\\rolek.brian\\Documents\\Projects\\CH2\\Data\\CoefficientPlots.csv",
         header=T)
omit <- cdat$Species[is.na(cdat$Covariate)]
cdat <- cdat[!cdat$Species %in% omit, ]
sp <- unique(cdat$Species)
pty <- c(0:6,8)
cols <- viridis(length(unique(cdat$Covariate)))
covs <- sort(table(cdat$Covariate), decreasing=T) 

fname<- "C:\\Users\\rolek.brian\\Documents\\Projects\\CH2\\docs\\coefplot.pdf"
pdf(fname, height=6, width=10) # width 36
par(mfrow=c(1,2), mar=c(3,3,0,0), oma=c(2, 10, 1, 0))
plot(NA, 
     ylim=c(length(sp)+0.25,0.75),
     xlim=c(-2.5, 1), 
     xaxt="n", yaxt="n")
axis(1, at= c(-2,-1, 0, 1))
axis(2, at=length(sp):1, labels=rev(sp), las=1)
mtext(side=1, "Coefficient value", line=3, cex=1.5)
mtext(side=2, "Species (common name)", line=11, cex=1.5)
abline(v=0, col="gray60", lwd=2, lty=2)
abline(h=0.5+0:length(sp), col="gray60", lwd=1, lty=1)
for (s in 1:length(sp)){
  cdat.temp <- cdat[cdat$Species==sp[s], ]
  cdat.temp <- cdat.temp[order(abs(cdat.temp$Mean), decreasing=T),]
  ylen <- dim(cdat.temp)[[1]]
  ymx <- ylen/2*0.1
  ycoord <- seq(-ymx, ymx, length=ylen)
  covnum <- match(cdat.temp$Covariate, rownames(covs))
  points(cdat.temp$Mean, s+ycoord, 
         pch=pty[covnum], col=cols[covnum], cex=2, lwd=2)
  for (i in 1:ylen){
    lines( cdat.temp[i,c(4,5)], c(s+ycoord[i],s+ycoord[i]),
           lwd=2, col=cols[covnum[i]])
  }
}
plot.new()
legend("left", 
       pch=c(pty,NA), col=c(cols, "black"), 
       legend=c(rownames(covs), "95% CIs"), title="Legend", 
       bty="n", pt.cex=2, pt.lwd=2, 
       lty=c(rep(NA, 8), 1), lwd=2,
       y.intersp=2, x.intersp=0.5)
dev.off()

#################
# Plot YSH
################
load("C:\\Users\\rolek.brian\\Documents\\Projects\\CH2\\Data\\DATA.RData")
load("E:\\chapter2\\outputs\\treat_GVS\\BBWA.RData")

w <- out$sims.list$w[,3]
beta.post <- out$sims.list$s.beta[w==1,8]
int.post <- out$sims.list$s.beta[w==1,1]
ysf <- datalfoc$tsh
x <- seq(min(ysf), max(ysf), length=100)
pred.x <- x*sd(stands.c$YrsSinceHarv, na.rm=T)+
          mean(stands.c$YrsSinceHarv, na.rm=T)
est <- array(NA, dim=c(length(beta.post), length(x)) )
for (i in 1:100){
est[,i] <- int.post + beta.post*x[i]
}
mns <- exp(apply(est, 2, mean))
lci <- exp(apply(est, 2, quantile, probs=0.025))
uci <- exp(apply(est, 2, quantile, probs=0.975))

# check for significance before plotting
mean(beta.post)
quantile(beta.post, prob=c(0.025, 0.975))

# plot
plot(pred.x, mns, type="n", 
     ylim=c(0,0.2), xlim=c(5,120),
     xaxt="n", yaxt="n")
polygon(x=c(pred.x, rev(pred.x)), y=c(lci, rev(uci)))
lines(pred.x, mns, lwd=3)
axis(1, at=c(20, 40, 60, 80, 100, 120), labels=c(20, NA, 60, NA, 100, NA))

