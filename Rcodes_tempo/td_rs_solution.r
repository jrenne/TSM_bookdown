# ==============================================================================
# Regime Switching Models
# ==============================================================================

#setwd("~/Google Drive/cours/ENSTA/R")

library(expm)


# ******************************************************************************
# Model:
#   X_t = mu_{s_t} + rho_{s_t}*X_{t-1} + sigma_{s_t} * epsilon_t
#   epsilon_t ~ N(0,1)
#   p_ij = P(s_t = j | s_t = i)
# ******************************************************************************


# ******************************************************************************
# Parameterization of Pi:

p_22 <- 0.95
p_11 <- 0.98
Pi <- matrix(c(p_11,1-p_22,1-p_11,p_22),2,2)


# ******************************************************************************
# Q1: Compute the undonditional distribution of s_t:
# ******************************************************************************


ergodic <- ((t(Pi))%^%1000)[,1]


# ******************************************************************************
# Parameterize mu, rho and Sigma:

rho2 <- 0.9
rho1 <- 0.7
rho <- matrix(c(rho1,rho2),1,2)
sigma2 <- .2
sigma1 <- .2
sigma <- matrix(c(sigma1,sigma2),1,2)
mu2 <- 0
mu1 <- .5
mu <- matrix(c(mu1,mu2),1,2)

unc_x <- (mu/(1-rho)) %*% ergodic

nb_sim <- 600
x <- matrix(0,nb_sim)
s <- matrix(0,nb_sim)



# ******************************************************************************
# Simulate s_t and X_t:

# Initialization
# =================
x[1] <- mu1/(1-rho1)
s[1] <- 1
x_coord <- 1
y_coord <- 1

for(t in 2:nb_sim){
  aux <- runif(1)
  if(s[t-1]==1){
    s[t] <- ifelse(aux<p_11,1,2)
  } else {
    s[t] <- ifelse(aux<p_22,2,1)
  }
  if(s[t]!=s[t-1]){
    x_coord <- c(x_coord,t,t)
    y_coord <- c(y_coord,s[t-1],s[t])
  }
  x[t] <- (s[t]==1)*(mu1 + rho1*x[t-1] + sigma1*rnorm(1)) +
    (s[t]==2)*(mu2 + rho2*x[t-1] + sigma2*rnorm(1))
}

# ******************************************************************************
# Compare histogram with the one of a Gaussian variable:
par(mfrow=c(1,1))
par(plt=c(0.1,0.9,0.1,0.9))
plot(density(x,bw=0.2),col="blue",lwd=2,main="Density of X(t) (in blue) compared with the Normal one (in red)")
x_norm <- mean(x)+sd(x)*rnorm(1000)
lines(density(x_norm,bw=0.2),col="red",lwd=2)


# ******************************************************************************
# Simulate trajectory:
par(mfrow=c(1,1))
par(plt=c(0.1,0.9,0.1,0.9))
plot(x,type='l',ylim=c(min(x),max(x)))
if(y_coord[length(y_coord)]==2){
  x_coord <- c(x_coord,nb_sim,nb_sim)
  y_coord <- c(y_coord,2,1)
} else {
  x_coord <- c(x_coord,nb_sim)
  y_coord <- c(y_coord,1)  
}
y_coord[y_coord==1] <- min(x)
y_coord[y_coord==2] <- max(x)
polygon(x_coord,y_coord,col='light grey',border = NA)
lines(x,lwd=2,type='l',ylim=c(min(x),max(x)),ljoin="round")


# ******************************************************************************
# Q2: Complete the following function, which returns the matrix eta,
#      whose lines are the eta_t, whose i-th entry is: f(x_t|x_{t-1},s_t = i)
# ******************************************************************************

compute_eta <- function(x,mu,sigma,rho,x_1){
  # N is the number of regimes:
  N <- max(dim(mu))
  T <- dim(x)[1]
  Ones <- matrix(1,T,1)
  y <- x%*%matrix(1,1,N) - Ones%*%mu - x_1%*%rho
  eta <- 1/sqrt(2*pi)*Ones%*%(1/sigma)*
    exp(-1/2*y^2/(Ones%*%(sigma^2)))
  return(eta)
}

x_1 <- as.matrix(x[1:(nb_sim-1)])
x <- as.matrix(x[2:nb_sim])
print(compute_eta(x,mu,sigma,rho,x_1))


# ******************************************************************************
# Q3: Complete the following function, which implements
#     the Kitagawa-Hamilton filter.
# ******************************************************************************

filter_prob <- function(eta,Pi,xi_0){
  T <- dim(eta)[1]
  N <- dim(eta)[2]
  xi <- matrix(0,T,N)
  xi[1,] <- (t(Pi)%*%xi_0)*eta[1,] /
    sum((t(Pi)%*%xi_0)*eta[1,])
  logl <- log(sum((t(Pi)%*%xi_0)*eta[1,]))
  for(t in 2:T){
    aux <- sum((t(Pi)%*%xi[t-1,])*eta[t,])
    xi[t,] <- (t(Pi)%*%xi[t-1,])*eta[t,] / aux
    logl <- logl + log(aux)
  }
  return(list(xi=xi,logl=logl))
}

x_1 <- x
x_1[2:NROW(x)] <- x[1:(NROW(x)-1),]
x_1[1] <- unc_x



# ******************************************************************************
# Q4: Use the filter to recover xi, the matrix of
#     filtered probabilities.
# ******************************************************************************

eta <- compute_eta(x,mu,sigma,rho,x_1)
#xi_0 <- uncond
xi_0 <- c(0,1)
xi   <- filter_prob(eta,Pi,xi_0)$xi

#stop()

par(mfrow=c(2,1))
par(plt=c(0.1,0.9,0.2,0.8))
plot(x,type='l',ylim=c(min(x),max(x)),main="X(t)",ylab="",xlab="")
if(y_coord[length(y_coord)]==2){
  x_coord <- c(x_coord,nb_sim,nb_sim)
  y_coord <- c(y_coord,2,1)
} else {
  x_coord <- c(x_coord,nb_sim)
  y_coord <- c(y_coord,1)  
}
y_coord[y_coord==1] <- min(x)
y_coord[y_coord==2] <- max(x)
polygon(x_coord,y_coord,col='light grey',border = NA)
lines(x,lwd=2,type='l',ylim=c(min(x),max(x)),ljoin="round")

plot(xi[,1],type='l',col="blue",lwd=2,main="Probability of being in the first regime (filtered: blue, smoothed:red)",ylab="",xlab="",ylim=c(0,1))



# ******************************************************************************
# Assume you do not know the model parameters.
# Q5: Complete the following function and maximize it to
#     get estimates of the model parameters
# ******************************************************************************

Loglik <- function(theta){
  mu <- matrix(theta[1:2],1,2)
  rho <- matrix(theta[3:4],1,2)
  sigma <- matrix(theta[5:6],1,2)
  p_22 <- exp(theta[7])/(1+exp(theta[7]))
  p_11 <- exp(theta[8])/(1+exp(theta[8]))
  Pi <- matrix(c(p_11,1-p_22,1-p_11,p_22),2,2)
  eta <- compute_eta(x,mu,sigma,rho,x_1)
  res <- filter_prob(eta,Pi,xi_0)
  return(- res$logl)
}

theta0 <- c(mu,rho,sigma,3,3)
#theta0 <- c(0*mu,0*rho,c(1,1))
print(Loglik(theta0))

opt1 <- optim(theta0,Loglik,
              method="Nelder-Mead",
              control=list(trace=TRUE,maxit=500),
              hessian=FALSE)
theta0 <- opt1$par
opt1 <- optim(theta0,Loglik,
              method="BFGS",
              control=list(trace=TRUE,maxit=100),
              hessian=FALSE)
theta0 <- opt1$par
opt1 <- optim(theta0,Loglik,
              method="Nelder-Mead",
              control=list(trace=TRUE,maxit=500),
              hessian=TRUE)
theta0 <- opt1$par

H <- opt1$hessian
H_1 <- solve(H)

A <- matrix(0,8,2)
A[,1] <- theta0
A[,2] <- sqrt(diag(H_1))
print(A)

