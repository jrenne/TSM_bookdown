
library(purrr)
library(expm)


par(plt=c(.1,.9,.1,.9))



simul.ARG <- function(nb.sim,mu,nu,rho,alpha=0,w0=NaN){
  unc.mean <- (alpha + nu) * mu/(1-rho)
  
  if(is.na(w0)){
    w_0 <- unc.mean
  }else{
    w_0 <- w0
  }
  W <- w_0
  w <- w_0
  for(t in 2:nb.sim){
    z <- rpois(1,rho*w/mu + alpha)
    w <- mu * rgamma(1,shape=nu+z)
    W <- c(W,w)
  }
  return(W)
}

W <- simul.ARG(300,mu=.5,nu=0,rho=.9,alpha=.1)

plot(W,type="l")


simul.compound.poisson <- function(nb.sim,Gamma,Pi,lambda,w0=NaN){
  
  if(is.na(w0)){
    w_0 <- 1 * Gamma
  }else{
    w_0 <- w0
  }
  
  W <- w_0
  w <- w_0
  for(t in 2:nb.sim){
    z <- sum(rbernoulli(n = w/Gamma, p = Pi))
    eps <- rpois(1,lambda)
    w <- Gamma * (z + eps)
    W <- c(W,w)
  }
  return(W)
}

nb.sim <- 200
Gamma <- 10
Pi <- .8
lambda <- .05
W <- simul.compound.poisson(nb.sim,Gamma,Pi,lambda)
plot(W,type="l")


simul.MS.AR <- function(nb.sim,mu.1,mu.2,rho.1,rho.2,sigma.1,sigma.2,P,w0=NaN){
  # s valued in {1,2}
  # p(1,2) = P( s_{t}=2 | s_{t-1}=1 )
  
  max.2.power <- 8
  Pk <- t(P)
  for(i in 1:max.2.power){
    Pk <- Pk %*% Pk
  }
  if(Pk[1,1]<.5){
    s <- 2
  }else{
    s <- 1
  }
  
  S <- s

  if(is.na(w0)){
    w_0 <- Pk[1,1] * mu.1/(1-rho.1) + Pk[2,1] * mu.2/(1-rho.2)
  }else{
    w_0 <- w0
  }
  
  W <- w_0
  w <- w_0
  
  for(t in 2:nb.sim){
    u <- runif(1)
    if(s==1){
      if(u>P[1,1]){
        s <- 2
      }
    }else{
      if(u<P[2,1]){
        s <- 1
      }
    }
    S <- c(S,s)
    
    if(s==1){
      w <- mu.1 + rho.1 * w + sigma.1 * rnorm(1)
    }else{
      w <- mu.2 + rho.2 * w + sigma.2 * rnorm(1)
    }
    
    W <- c(W,w)
  }
  return(list(S=S,W=W,Pk=Pk))
}

P <- matrix(0,2,2)
P[1,] <- c(.95,.05)
P[2,] <- c(.01,.99)

nb.sim <- 500
mu.1 <- 1
mu.2 <- 0
rho.1 <- .9
rho.2 <- .9
sigma.1 <- .1
sigma.2 <- .02

WS <- simul.MS.AR(nb.sim,mu.1,mu.2,
                  rho.1,rho.2,
                  sigma.1,sigma.2,
                  P,w0=NaN)

par(mfrow=c(2,1))
plot(WS$W,type="l")
plot(WS$S,type="l")
