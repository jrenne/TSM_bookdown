
setwd("~/Google Drive/cours/UNIL/")
source("Rcode/various_proc_TS.R")
setwd("~/Google Drive/EDF_partenariat/Slides/")

r.bar <- .03
sigma.bar <- .012

phi <- .98
mu <- r.bar * (1-phi)
sigma2 <- 16 * sigma.bar^2 * (1 - phi^2)
sigma2 <- sigma.bar^2 * (1 - phi^2)

X.bar <- mu/(1-phi)


Model <- list()

Model$phi <- phi
Model$mu <- mu
Model$Sigma <- sqrt(sigma2)/12
Model$delta.0 <- 0
Model$delta.1 <- 1/12

compute.ab <- function(Model,h){
  phi <- Model$phi
  mu <- Model$mu
  Sigma <- Model$Sigma
  delta.0 <- Model$delta.0
  delta.1 <- Model$delta.1
  
  A.n.bar <- - delta.0
  B.n.bar <- - delta.1
  all.A <- - A.n.bar
  all.B <- - B.n.bar
  for(n in 2:h){
    A.n.bar <- -delta.0 + A.n.bar + t(B.n.bar) %*% mu +
      .5 * t(B.n.bar) %*% Sigma %*% t(Sigma) %*% B.n.bar
    B.n.bar <- t(- t(delta.1) + t(B.n.bar) %*% phi)
    all.A <- c(all.A,- A.n.bar/n)
    all.B <- cbind(all.B,- B.n.bar/n)
  }
  return(list(
    all.A = all.A,
    all.B = all.B
  ))
}




stop()





h <- 120

loadings <- compute.ab(Model,h)

all.A.H <- loadings$all.A
all.B.H <- loadings$all.B

par(mfrow=c(1,1))
plot(12*t(all.A.H + all.B.H * X.bar),ylim=c(.01,.05),type="l")



Model.star <- Model
Model.star$phi <- .99
Model.star$mu <- .05 * (1 - Model.star$phi)

loadings <- compute.ab(Model.star,h)

all.A <- loadings$all.A
all.B <- loadings$all.B

# Unconditional TS:
par(mfrow=c(1,1))
plot(t(all.A + all.B * X.bar),ylim=c(.02,.05),type="l")

# biased Model:
Model.biased <- Model
Model.biased$phi <- .9
Model.biased$mu <- X.bar * (1 - Model.biased$phi)

loadings <- compute.ab(Model.biased,h)

all.A.biased <- loadings$all.A
all.B.biased <- loadings$all.B


T <- 200
nb.sim <- 1
x.sim <- sim.arma(mu,phi,theta=c(1),sqrt(sigma2),T,X.bar,nb.sim)
i.t <- 12*(Model$delta.0 + Model$delta.1 * x.sim)

matur <- 120
i.h <- 12*( all.A[matur]+all.B[matur]*x.sim)
i.h.H <- 12*(all.A.H[matur]+all.B.H[matur]*x.sim)
i.h.H.biased <- 12*(all.A.biased[matur]+all.B.biased[matur]*x.sim)

FILE = "/figures/Figure_atsm_downward.pdf"
pdf(file=paste(getwd(),FILE,sep=""),pointsize=10,width=10, height=7)
par(mfrow=c(2,1))
par(plt=c(.15,.95,.15,.95))

plot(i.t,type="l",lwd=2,
     ylim=c(min(i.t,i.h,i.h.H),max(i.t,i.h,i.h.H)+.02),
     xlab="",ylab="")
lines(i.h, col="red",lwd=2)
lines(i.h.H, col="blue",lwd=2)
lines(i.h.H.biased, col="dark grey",lwd=2,lty=2)
abline(h=X.bar,col="grey")

legend("topright", # places a legend at the appropriate place c("Health","Defense"), # puts text in the legend 
       c(
         "Short-term rate",
         "10-year rate",
         "10-year rate if no term-premium",
         "10-year rate if no term-premium, biased"
       ),
       #pt.bg = c("black","red","blue"),
       #lty=c(1,1), # gives the legend appropriate symbols (lines)       
       lwd=c(2), # line width
       lty=c(1,1,1,2),
       #pch=c(19,19,19),
       col=c("black","red","blue","dark grey"), # gives the legend lines the correct color and width
       seg.len = 4,
       bg = "white"
)


plot(i.h - i.h.H,col="blue",type="l",
     ylim=c(min(i.h - i.h.H,i.h - i.h.H.biased),max(i.h - i.h.H,i.h - i.h.H.biased)),lwd=2,
     xlab="",ylab="")
lines(i.h - i.h.H.biased,col="dark grey",lwd=2,lty=2)

legend("topright", # places a legend at the appropriate place c("Health","Defense"), # puts text in the legend 
       c(
         "10-year term-premium ('True')",
         "10-year term-premium, biased"
       ),
       #pt.bg = c("black","red","blue"),
       #lty=c(1,1), # gives the legend appropriate symbols (lines)       
       lwd=c(2), # line width
       lty=c(1,2),
       #pch=c(19,19,19),
       col=c("blue","dark grey"), # gives the legend lines the correct color and width
       seg.len = 4,
       bg = "white"
)


dev.off()

stop()




# First attempt with CCAPM:

phi <- 0.5
mu <- 0.01 * (1 - phi)
sigma.bar <- .01
sigma2 <- sigma.bar^2 * (1 - phi^2)
sigma.c <- sqrt(sigma2)
  
gamma <- 7.5
beta <- 0.999

lambda.0 <- gamma * sigma.c

mu.star <- mu - sigma.c * lambda.0

Model$phi <- phi
Model$mu <- mu.star
Model$Sigma <- sigma.c
Model$delta.0 <- -log(beta) + mu * rho *(1 - phi) - rho^2 * sigma.c^2 / 2
Model$delta.1 <- rho * phi

loadings <- compute.ab(Model,h)

all.A <- loadings$all.A
all.B <- loadings$all.B

r.bar <- Model$delta.0 + Model$delta.1 * mu/(1-phi)
  
par(mfrow=c(1,1))
plot(t(all.A + all.B * r.bar),type="l")

stop()


# Multivariate case:

Model$phi <- matrix(c(.9,.1,0,.9),2,2)
Model$mu <- c(.001,.002)
Model$Sigma <- diag(rep(sqrt(sigma2),2))
Model$delta.0 <- 0
Model$delta.1 <- matrix(c(1,1),2,1)

loadings <- compute.ab(Model,10)



