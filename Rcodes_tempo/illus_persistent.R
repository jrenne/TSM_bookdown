# =========================================
# Code for the nonstationarity lesson
# =========================================

rho <- 1


# Clear environment and console
rm(list = ls(all = TRUE)) # clear environment 

setwd("~/Google Drive/cours/UNIL/")
source("Rcode/various_proc_TS.R")
setwd("~/Google Drive/EDF_partenariat/Slides/")

library(chron)


# ================================
# Examples of MA processes


y.0 <- c(0)
c <- 0
nb.sim <- 100
theta <- c(1)
phi <- c(1)
sigma <- 1


vec.T <- c(50,400)

FILE = "/figures/Figure_NonStat_Probl1.pdf"
pdf(file=paste(getwd(),FILE,sep=""),pointsize=10,width=8, height=4)
par(mfrow=c(1,2))
par(plt=c(.2,.9,.2,.85))

for(T in vec.T){
  y.sim <- sim.arma(c,phi,theta,sigma,T,y.0,nb.sim)
  
  vec.rho.hat <- NULL
  for(i in 1:nb.sim){
    y.t <- y.sim[2:T,i]
    y.t_1 <- y.sim[1:(T-1),i]
    eq <- lm(y.t ~ y.t_1)
    rho.hat <- eq$coefficients[2]
    vec.rho.hat <- c(vec.rho.hat,rho.hat)
  }
  dens <- density(vec.rho.hat)
  
  main.t <- paste("T =",toString(T))
  
  mean.rho.hat <- sum(dens$x*dens$y)/sum(dens$y)
  plot(dens,xlab="Estimated phi",ylab="",
       main=main.t)
  abline(v=mean.rho.hat,col="red")
  abline(v=1-5.3/T,col="blue")
  
}

dev.off()




vec.T <- c(50,400)

FILE = "/figures/Figure_NonStat_Probl2.pdf"
pdf(file=paste(getwd(),FILE,sep=""),pointsize=10,width=10, height=4)
par(mfrow=c(1,2))
par(plt=c(.2,.9,.2,.85))

for(T in vec.T){
  
  y.sim <- sim.arma(c,phi,theta,sigma,T,y.0,nb.sim)
  x.sim <- sim.arma(c,phi,theta,sigma,T,y.0,nb.sim)
  
  #y.sim <- 2*x.sim + 3*matrix(rnorm(T*nb.sim),T,nb.sim)
  
  vec.rho.hat <- NULL
  vec.stdv <- NULL
  for(i in 1:nb.sim){
    y.t <- y.sim[,i]
    x.t <- x.sim[,i]
    eq <- lm(y.t ~ x.t)
    rho.hat <- eq$coefficients[2]
    vec.rho.hat <- c(vec.rho.hat,rho.hat)
    vec.stdv <- c(vec.stdv,sqrt(vcov(eq)[2,2]))
  }
  dens <- density(vec.rho.hat)
  
  main.t <- paste("T=",toString(T),", Mean of OLS-based std dev of phi:",toString(round(mean(vec.stdv),3)),sep="")
  
  mean.rho.hat <- sum(dens$x*dens$y)/sum(dens$y)
  plot(dens,xlab="Estimated phi",ylab="",lwd=2,main=main.t)
  abline(v=0,col="red")
  #abline(v=+2*mean(vec.stdv))
  #lines(seq(-2,2,by=.001),dnorm(seq(-2,2,by=.001),mean = 0,sd = mean(vec.stdv)),col="red")
  
}

dev.off()


nb.sim <- 1000

k <- 12

vec.T <- c(50,400)

FILE = "/figures/Figure_NonStat_Probl3.pdf"
pdf(file=paste(getwd(),FILE,sep=""),pointsize=10,width=9, height=6)
par(mfrow=c(2,3))
par(plt=c(.2,.9,.2,.85))

phi.y <- 0
phi.x <- .9

for(T in vec.T){
  
  y.sim0 <- sim.arma(c,phi.y,theta,sigma,T+(k-1),y.0,nb.sim)
  y.sim <- y.sim0[1:T,]
  for(i in 2:k){
    y.sim <- y.sim + y.sim0[i:(T+i-1),]
  }

  x.sim <- sim.arma(c,phi.x,theta,sigma,T,y.0,nb.sim)
  
  #y.sim <- 2*x.sim + 3*matrix(rnorm(T*nb.sim),T,nb.sim)
  
  vec.rho.hat <- NULL
  vec.stdv <- NULL
  vec.tstud <- NULL
  vec.R2 <- NULL
  for(i in 1:nb.sim){
    y.t <- y.sim[,i]
    x.t <- x.sim[,i]
    eq <- lm(y.t ~ x.t)
    rho.hat <- eq$coefficients[2]
    vec.rho.hat <- c(vec.rho.hat,rho.hat)
    vec.stdv <- c(vec.stdv,sqrt(vcov(eq)[2,2]))
    vec.tstud <- c(vec.tstud,rho.hat/sqrt(vcov(eq)[2,2]))
    vec.R2 <- c(vec.R2,summary(eq)$r.sq)
  }
  dens <- density(vec.rho.hat)
  
  main.t <- paste("T=",toString(T),", Mean of OLS-based std dev of phi:",toString(round(mean(vec.stdv),3)),sep="")
  mean.rho.hat <- sum(dens$x*dens$y)/sum(dens$y)
  plot(dens,xlab=expression(phi),
       ylab="",lwd=2,main=main.t)
  abline(v=0,col="red")
  
  dens <- density(vec.R2)
  main.t <- paste("T=",toString(T),", R2",sep="")
  mean.rho.hat <- sum(dens$x*dens$y)/sum(dens$y)
  plot(dens,xlab="R2",ylab="",lwd=2,main=main.t)
  
  dens <- density(vec.tstud)
  main.t <- paste("T=",toString(T),", Student t",sep="")
  mean.rho.hat <- sum(dens$x*dens$y)/sum(dens$y)
  plot(dens,xlab="Student t",ylab="",lwd=2,main=main.t)
  #abline(v=+2*mean(vec.stdv))
  #lines(seq(-2,2,by=.001),dnorm(seq(-2,2,by=.001),mean = 0,sd = mean(vec.stdv)),col="red")
  
}

dev.off()



