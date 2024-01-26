
setwd("~/Google Drive/EDF_partenariat/Slides")


# Illustration de RRA

C.l <- .8
C.h <- 1.2
cc <- c(C.l,C.h)

FILE = "/figures/Figure_RRA_IES.pdf"
#pdf(file=paste(getwd(),FILE,sep=""),pointsize=10,width=6, height=5)  
pdf(file=paste(getwd(),FILE,sep=""),pointsize=10,width=7, height=5)
par(mfrow=c(1,2))
par(plt=c(.2,.9,.2,.85))

C <- seq(.7,1.5,by=.01)

gamma <- 2
E.u <- mean(cc^(1-gamma)/(1-gamma))
c.CE <- (E.u*(1-gamma))^(1/(1-gamma))
U <- C^(1-gamma)/(1-gamma)
plot(C,U,type="l",lwd=2,
     xlab="Consumption",ylab="Utility",
     main=expression(paste(gamma,"=2",sep="")))
points(cc,cc^(1-gamma)/(1-gamma),pch=19,col="blue",cex=1.5)
lines(cc,cc^(1-gamma)/(1-gamma),col="blue")
lines(c(c.CE,mean(cc)),c(E.u,E.u),lty=1,col="dark grey",lwd=1)
lines(c(c.CE,c.CE),c(-1000,E.u),lty=1,col="dark grey",lwd=1)
lines(c(mean(cc),mean(cc)),c(-1000,1/(1-gamma)),lty=1,col="dark grey",lwd=1)
points(c(1,1),c(mean(cc^(1-gamma)/(1-gamma)),1/(1-gamma)),col="red",pch=19)
lines(c(1,1),c(mean(cc^(1-gamma)/(1-gamma)),1/(1-gamma)),col="red")
rp <- - mean(cc^(1-gamma)/(1-gamma)) + 1/(1-gamma)
text(x = 1.2,y=1/(1-gamma) - .2,
     labels = paste("Risk premium = ",toString(round(mean(cc)-c.CE,2)),sep=""),
     col="red")
lines(c(mean(c(c.CE,1)),1.2),c(C[1]^(1-gamma)/(1-gamma),1/(1-gamma) - .22),col="red")

gamma <- 10
E.u <- mean(cc^(1-gamma)/(1-gamma))
c.CE <- (E.u*(1-gamma))^(1/(1-gamma))
U <- C^(1-gamma)/(1-gamma)
plot(C,U,type="l",lwd=2,
     xlab="Consumption",ylab="Utility",
     main=expression(paste(gamma,"=10",sep="")))
points(cc,cc^(1-gamma)/(1-gamma),pch=19,col="blue",cex=1.5)
lines(cc,cc^(1-gamma)/(1-gamma),col="blue")
lines(c(c.CE,mean(cc)),c(E.u,E.u),lty=1,col="dark grey",lwd=1)
lines(c(c.CE,c.CE),c(-1000,E.u),lty=1,col="dark grey",lwd=1)
lines(c(mean(cc),mean(cc)),c(-1000,1/(1-gamma)),lty=1,col="dark grey",lwd=1)
points(c(1,1),c(mean(cc^(1-gamma)/(1-gamma)),1/(1-gamma)),col="red",pch=19)
lines(c(1,1),c(mean(cc^(1-gamma)/(1-gamma)),1/(1-gamma)),col="red")
rp <- - mean(cc^(1-gamma)/(1-gamma)) + 1/(1-gamma)
text(x = 1.2,y=1/(1-gamma) - 1,
     labels = paste("Risk premium = ",toString(round(mean(cc)-c.CE,2)),sep=""),
     col="red")
lines(c(mean(c(c.CE,1)),1.2),c(C[1]^(1-gamma)/(1-gamma),1/(1-gamma) - 1.05),col="red")

dev.off()


# Illustration de influence news.

beta <- .9
gamma <- seq(0.5001,10,by=.1)

FILE = "/figures/Figure_EZ_news.pdf"
#pdf(file=paste(getwd(),FILE,sep=""),pointsize=10,width=6, height=5)  
pdf(file=paste(getwd(),FILE,sep=""),pointsize=10,width=10, height=6)
par(mfrow=c(1,2))
par(plt=c(.15,.9,.2,.8))

omega <- .05
E.m.0 <- - log(0.5 * exp(beta*(1-gamma)*omega) + 0.5 * exp(-beta*(1-gamma)*omega))
P <- .5 * exp(E.m.0 - beta*(1-gamma)*omega)
P.check <- .5 * exp(E.m.0 + beta*(1-gamma)*omega) + .5 * exp(E.m.0 - beta*(1-gamma)*omega) # should be equal to 1
plot(gamma,P,type="l",lwd=2,
     xlab=expression(gamma),ylab="",
     main=expression(paste("Price of an asset that pays 1 under Case II (",omega,"=0.05)",sep="")))
abline(v=1,col="red")
abline(h=0.5,col="blue",lty=2)

omega <- .2
E.m.0 <- - log(0.5 * exp(beta*(1-gamma)*omega) + 0.5 * exp(-beta*(1-gamma)*omega))
P <- .5 * exp(E.m.0 - beta*(1-gamma)*omega)
P.check <- .5 * exp(E.m.0 + beta*(1-gamma)*omega) + .5 * exp(E.m.0 - beta*(1-gamma)*omega) # should be equal to 1
plot(gamma,P,type="l",lwd=2,
     xlab=expression(gamma),ylab="",
     main=expression(paste("Price of an asset that pays 1 under Case II (",omega,"=0.2)",sep="")))
abline(v=1,col="red")
abline(h=0.5,col="blue",lty=2)

dev.off()




# Illustration of early resolution of uncertainty.

FILE = "/figures/Figure_EZ_early.pdf"
#pdf(file=paste(getwd(),FILE,sep=""),pointsize=10,width=6, height=5)  
pdf(file=paste(getwd(),FILE,sep=""),pointsize=10,width=8, height=4)
par(mfrow=c(1,3))
par(plt=c(.15,.9,.2,.8))

for(rho in c(.5,.9,1.5)){

all.gamma.values <- seq(.1001,6,by=.1)
all.C.minus.D.values <- NULL

for(gamma in all.gamma.values){
  c.h <- 1.2
  c.l <- .8
  beta <- .5
  
  U.2.h <- (1 - beta)^(1/(1-rho)) * c.h
  U.2.l <- (1 - beta)^(1/(1-rho)) * c.l
  
  U.1.C.hh <- ((1-beta)*c.h^(1-rho) + beta * U.2.h^(1-rho))^(1/(1-rho))
  U.1.C.hl <- ((1-beta)*c.h^(1-rho) + beta * U.2.l^(1-rho))^(1/(1-rho))
  U.1.C.lh <- ((1-beta)*c.l^(1-rho) + beta * U.2.h^(1-rho))^(1/(1-rho))
  U.1.C.ll <- ((1-beta)*c.l^(1-rho) + beta * U.2.l^(1-rho))^(1/(1-rho))
  
  U.0.C <- (1 - beta)^(1/(1-rho)) * (.25 * U.1.C.hh^(1-gamma) +.25 * U.1.C.hl^(1-gamma) +
                                       .25 * U.1.C.lh^(1-gamma) +.25 * U.1.C.ll^(1-gamma)
  )^(1/(1-gamma))
  
  E.U2gamma <- .5 * U.2.h^(1 - gamma) + .5 * U.2.l^(1 - gamma)
  
  U.1.D.h <- ((1-beta)*c.h^(1-rho) + beta * E.U2gamma^((1-rho)/(1-gamma)))^(1/(1-rho))
  U.1.D.l <- ((1-beta)*c.l^(1-rho) + beta * E.U2gamma^((1-rho)/(1-gamma)))^(1/(1-rho))
  
  U.0.D <- (1 - beta)^(1/(1-rho)) * (.5 * U.1.D.h^(1-gamma) +.5 * U.1.D.l^(1-gamma))^(1/(1-gamma))
  
  all.C.minus.D.values <- c(all.C.minus.D.values,U.0.C-U.0.D)
}

if(rho==.5){
  plot(all.gamma.values,all.C.minus.D.values,type="l",lwd=2,
       xlab=expression(gamma),ylab="",
       main=expression(paste("Difference between ",U[0]^C," and ",U[0]^D," when ",rho,"=0.50",sep=""))
  )
}
if(rho==.9){
  plot(all.gamma.values,all.C.minus.D.values,type="l",lwd=2,
       xlab=expression(gamma),ylab="",
       main=expression(paste("Difference between ",U[0]^C," and ",U[0]^D," when ",rho,"=0.90",sep=""))
  )
}
if(rho==1.5){
  plot(all.gamma.values,all.C.minus.D.values,type="l",lwd=2,
       xlab=expression(gamma),ylab="",
       main=expression(paste("Difference between ",U[0]^C," and ",U[0]^D," when ",rho,"=1.50",sep=""))
  )
}
abline(v=rho,col="red")
abline(h=0,col="blue",lty=2)
}
dev.off()



FILE = "/figures/Figure_PowerU_smoothing.pdf"
#pdf(file=paste(getwd(),FILE,sep=""),pointsize=10,width=6, height=5)  
pdf(file=paste(getwd(),FILE,sep=""),pointsize=10,width=6, height=3)
par(mfrow=c(1,1))
par(plt=c(.3,.7,.3,.95))

# Illustration of IES and smoothing
W <- 1
delta <- 1
R <- .05

gamma <- 1.5
C1 <- W * (1 + R) / ((1/delta * 1/(1+R))^(-1/gamma) + (1+R))
C2 <- (W - C1) * (1 + R)
plot(c(C1,C2),pch=19,xlab="Period",ylab="Consumption",ylim=c(.50,.53))
lines(c(C1,C2))
gamma <- 5
C1 <- W * (1 + R) / ((1/delta * 1/(1+R))^(-1/gamma) + (1+R))
C2 <- (W - C1) * (1 + R)
points(c(C1,C2),pch=19,col="red")
lines(c(C1,C2),col="red")
gamma <- 50
C1 <- W * (1 + R) / ((1/delta * 1/(1+R))^(-1/gamma) + (1+R))
C2 <- (W - C1) * (1 + R)
points(c(C1,C2),pch=19,col="blue")
lines(c(C1,C2),col="blue")

legend("topright", # places a legend at the appropriate place c("Health","Defense"), # puts text in the legend 
       c(expression(paste(gamma,"=1.5"," (IES=2/3)",sep="")),
         expression(paste(gamma,"=5"," (IES=1/5)",sep="")),
         expression(paste(gamma,"=50"," (IES=1/50)",sep=""))
       ),
       pt.bg = c("black","red","blue"),
       lty=c(1,1), # gives the legend appropriate symbols (lines)       
       lwd=c(2,2), # line width
       pch=c(19,19,19),
       col=c("black","red","blue"), # gives the legend lines the correct color and width
       seg.len = 4
)

dev.off()





# Illustration of habit

delta <- 1
lambda <- .49
rho <- .5
gamma <- 3

r.f <- c(rep(0.07,5),c(-.2,-.2,-.2))

T <- length(r.f)
C <- 1
C.t <- C
X <- lambda * C /(1-rho)
#X <- lambda
X.t <- X
U <- delta * (C.t - X.t)^(1 - gamma)/(1 - gamma)
U.no.habit <- delta * (C.t)^(1 - gamma)/(1 - gamma)
for(t in 1:(T)){
  C.t <- 1/(1-lambda) * ((1/delta * 1/(1+r.f[t]))^(-1/gamma) * (C.t - X.t) + rho * X.t)
  C <- c(C,C.t)
  X.t <- rho * X.t + lambda * C.t
  X <- c(X,X.t)
  U <- U + delta * (C.t - X.t)^(1 - gamma)/(1 - gamma)
  U.no.habit <- U.no.habit + delta * (C.t)^(1 - gamma)/(1 - gamma)
  #print(c(t, delta * (C.t)^(1 - gamma)/(1 - gamma)))
}
U1 <- U
U1.no.habit <- U.no.habit
C1 <- C
X1 <- X
R.f.1 <- r.f

par(mfrow=c(2,3))
plot(r.f,type="l")
plot(C,type="l",ylim=c(min(X),max(C)))
lines(X,col="blue")
plot(C-X,type="l")

r.f <- c(rep(0.025,8))

T <- length(r.f)
C <- 1
C.t <- C
X <- lambda * C /(1-rho)
#X <- lambda
X.t <- X
U <- delta * (C.t - X.t)^(1 - gamma)/(1 - gamma)
U.no.habit <- delta * (C.t)^(1 - gamma)/(1 - gamma)
for(t in 1:(T)){
  C.t <- 1/(1-lambda) * ((1/delta * 1/(1+r.f[t]))^(-1/gamma) * (C.t - X.t) + rho * X.t)
  C <- c(C,C.t)
  X.t <- rho * X.t + lambda * C.t
  X <- c(X,X.t)
  U <- U + delta * (C.t - X.t)^(1 - gamma)/(1 - gamma)
  U.no.habit <- U.no.habit + delta * (C.t)^(1 - gamma)/(1 - gamma)
  #print(c(t, delta * (C.t)^(1 - gamma)/(1 - gamma)))
}
U2 <- U
U2.no.habit <- U.no.habit
C2 <- C
X2 <- X
R.f.2 <- r.f

plot(r.f,type="l")
plot(C2,type="l",ylim=c(min(X),max(C1,C2)))
lines(C1,col="red")
lines(X,col="blue")
plot(C2-X2,type="l",ylim=c(min(C1-X1,C2-X2),max(C1-X1,C2-X2)))
lines(C1-X1,col="red")

print(c(U1,U2))
print(c(U1.no.habit,U2.no.habit))



FILE = "/figures/Figure_habit.pdf"
#pdf(file=paste(getwd(),FILE,sep=""),pointsize=10,width=6, height=5)  
pdf(file=paste(getwd(),FILE,sep=""),pointsize=10,width=7, height=3)
par(mfrow=c(1,3))
par(plt=c(.15,.9,.15,.85))

plot(R.f.1,col="blue",lwd=2,pch=19,ylim=c(-.25,.13),
     main=expression(paste("Short-term interest rate ",R[f*","*t],sep="")),xlab="time",ylab="")
lines(R.f.1,col="blue")
points(R.f.2,col="red",pch=19)
lines(R.f.2,col="red")


legend("topright", # places a legend at the appropriate place c("Health","Defense"), # puts text in the legend 
       c("Scenario A","Scenario B"),
       lty=c(1,1), # gives the legend appropriate symbols (lines)       
       lwd=c(2,2), # line width
       pch=c(19,19),
       col=c("blue","red"), # gives the legend lines the correct color and width
       seg.len = 4
)

plot(C1,col="blue",lwd=2,pch=19,ylim=c(1,1.01),
     main=expression(paste("Consumption ",C[t],sep="")),xlab="time",ylab="")
lines(C1,col="blue")
points(C2,col="red",pch=19)
lines(C2,col="red")

plot(C1-X1,col="blue",lwd=2,pch=19,ylim=c(0.01,0.025),
     main=expression(C[t]-X[t]),xlab="time",ylab="")
lines(C1-X1,col="blue")
points(C2-X2,col="red",pch=19)
lines(C2-X2,col="red")


dev.off()




