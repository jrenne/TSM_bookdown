#
# Illustration Kalman
#
# (only non-pricing things. The pricing things are in Shiny/run.kalman.R)


T <- 100
alpha1 <- .5
alpha2 <- .95
Alpha <- diag(c(alpha1,alpha2))
d_11 <- 1
d_12 <- .5
d_21 <- .5
d_22 <- 2

D <- matrix(c(d_11,d_21,d_12,d_22),2,2)
gamma1 <- 1
gamma2 <- 2
Gamma <- matrix(c(gamma1,gamma2),2,1)

phi <- .8

Y <- NULL
X <- NULL
Alpha.Y_1 <- NULL
y <- c(0,0)
x <- 0
for(i in 1:T){
  Alpha.Y_1 <- rbind(Alpha.Y_1,c(Alpha %*% y))
  y <- Alpha %*% y + Gamma * x + D %*% rnorm(2)
  x <- phi * x + rnorm(1)
  Y <- rbind(Y,t(y))
  X <- rbind(X,x)
}

Rfunction <- function(M,RHO,t=0){
 return(M %*% t(M))
}
Qfunction <- function(N,RHO,t=0){
 return(N %*% t(N))
}


nu_t <- matrix(0,T,1)
H <- phi
N <- 1
mu_t <- Alpha.Y_1
G <- Gamma
M <- D
Sigma_0 <- 1/(1-phi^2)
rho_0 <- 0
  
filter.res <-     Kalman_filter(Y,nu_t,H,N,mu_t,G,M,Sigma_0,rho_0,indic_pos=0)
smoother.res <- Kalman_smoother(Y,nu_t,H,N,mu_t,G,M,Sigma_0,rho_0,indic_pos=0)


# FILE = "/figures/Figure_Kalman1.pdf"
# pdf(file=paste(getwd(),FILE,sep=""),pointsize=10,width=6, height=6)
par(mfrow=c(3,1))
par(plt=c(.1,.95,.15,.85))

par(mfrow=c(3,1))
plot(Y[,1],type="l",main=expression(Y[1*","*t]),lwd=2,xlab="",ylab="")
plot(Y[,2],type="l",main=expression(Y[2*","*t]),lwd=2,xlab="",ylab="")
plot(X,type="l",main=expression(X[t]),lwd=2,xlab="",ylab="")

# legend("topright", # places a legend at the appropriate place c("Health","Defense"), # puts text in the legend 
#        c(
#          "Short-term rate",
#          "Term structure of real yields",
#          "Term structure of real yields (without risk premiums)"),
#        #pt.bg = c("black","red","blue"),
#        #lty=c(1,1), # gives the legend appropriate symbols (lines)       
#        lwd=c(2), # line width
#        lty=c(1,1,1),
#        #pch=c(19,19,19),
#        col=c("black","red","blue"), # gives the legend lines the correct color and width
#        seg.len = 4,
#        bg = "white"
# )

#dev.off()


# FILE = "/figures/Figure_Kalman2.pdf"
# pdf(file=paste(getwd(),FILE,sep=""),pointsize=10,width=8, height=6)
par(mfrow=c(2,1))
par(plt=c(.1,.95,.15,.8))

plot(X,type="l",main=expression(X[t]),lwd=2,xlab="",ylab="")
lines(filter.res$r,col="red",lwd=2)
lines(smoother.res$r,col="blue",lwd=2,lty=1)

legend("topright", # places a legend at the appropriate place c("Health","Defense"), # puts text in the legend 
       c(
         expression(X[t]),
         expression(paste("Filtered values of ",X[t],sep="")),
         expression(paste("Smoothed values of ",X[t],sep=""))),
       #pt.bg = c("black","red","blue"),
       #lty=c(1,1), # gives the legend appropriate symbols (lines)       
       lwd=c(2), # line width
       lty=c(1,1,1),
       #pch=c(19,19,19),
       col=c("black","red","blue"), # gives the legend lines the correct color and width
       seg.len = 4,
       bg = "white"
)

plot(sqrt(filter.res$Sigma_tt),type="l",main=expression(X[t]),lwd=2,xlab="",ylab="",col="red",
     ylim=c(min(sqrt(smoother.res$S_smooth)),max(sqrt(filter.res$Sigma_tt))))
lines(sqrt(smoother.res$S_smooth),lwd=2,col="blue")

legend("topright", # places a legend at the appropriate place c("Health","Defense"), # puts text in the legend 
       c(
         expression(paste("Filtering standard error ",sqrt(Var(X[t*"|"*t])),sep="")),
         expression(paste("Smoothing standard error ",sqrt(Var(X[t*"|"*T])),sep=""))),
       #pt.bg = c("black","red","blue"),
       #lty=c(1,1), # gives the legend appropriate symbols (lines)       
       lwd=c(2), # line width
       lty=c(1,1,3),
       #pch=c(19,19,19),
       col=c("black","red","red"), # gives the legend lines the correct color and width
       seg.len = 4,
       bg = "white"
)

# dev.off()


# Missing data issue: remove observations.

Y.modif <- Y
Y.modif[30:50,1] <- NaN
Y.modif[40:70,2] <- NaN

filter.res <-     Kalman_filter(Y.modif,nu_t,H,N,mu_t,G,M,Sigma_0,rho_0,indic_pos=0)
smoother.res <- Kalman_smoother(Y.modif,nu_t,H,N,mu_t,G,M,Sigma_0,rho_0,indic_pos=0)




# FILE = "/figures/Figure_Kalman3.pdf"
# pdf(file=paste(getwd(),FILE,sep=""),pointsize=10,width=6, height=6)
par(mfrow=c(3,1))
par(plt=c(.1,.95,.15,.85))

par(mfrow=c(3,1))
plot(Y.modif[,1],type="l",main=expression(Y[1*","*t]),lwd=2,xlab="",ylab="")
plot(Y.modif[,2],type="l",main=expression(Y[2*","*t]),lwd=2,xlab="",ylab="")
plot(X,type="l",main=expression(X[t]),lwd=2,xlab="",ylab="")

# legend("topright", # places a legend at the appropriate place c("Health","Defense"), # puts text in the legend 
#        c(
#          "Short-term rate",
#          "Term structure of real yields",
#          "Term structure of real yields (without risk premiums)"),
#        #pt.bg = c("black","red","blue"),
#        #lty=c(1,1), # gives the legend appropriate symbols (lines)       
#        lwd=c(2), # line width
#        lty=c(1,1,1),
#        #pch=c(19,19,19),
#        col=c("black","red","blue"), # gives the legend lines the correct color and width
#        seg.len = 4,
#        bg = "white"
# )

# dev.off()





# FILE = "/figures/Figure_Kalman3_new.pdf"
# pdf(file=paste(getwd(),FILE,sep=""),pointsize=10,width=4, height=2)
par(mfrow=c(3,1))
par(plt=c(.15,.95,.3,.95))

par(mfrow=c(3,1))
plot(Y.modif[,1],type="l",main=expression(y[1*","*t]),lwd=2,xlab="",ylab=expression(y[1*","*t]))
plot(Y.modif[,2],type="l",main=expression(y[2*","*t]),lwd=2,xlab="",ylab=expression(y[2*","*t]))
plot(X,type="l",main=expression(w[t]),lwd=2,xlab="",ylab=expression(w[t]))

# dev.off()



# FILE = "/figures/Figure_Kalman4.pdf"
# pdf(file=paste(getwd(),FILE,sep=""),pointsize=10,width=8, height=6)
par(mfrow=c(2,1))
par(plt=c(.1,.95,.15,.8))

plot(X,type="l",main=expression(X[t]),lwd=2,xlab="",ylab="")
lines(filter.res$r,col="red",lwd=2)
lines(smoother.res$r,col="blue",lwd=2,lty=1)

legend("topright", # places a legend at the appropriate place c("Health","Defense"), # puts text in the legend 
       c(
         expression(X[t]),
         expression(paste("Filtered values of ",X[t],sep="")),
         expression(paste("Smoothed values of ",X[t],sep=""))),
       #pt.bg = c("black","red","blue"),
       #lty=c(1,1), # gives the legend appropriate symbols (lines)       
       lwd=c(2), # line width
       lty=c(1,1,1),
       #pch=c(19,19,19),
       col=c("black","red","blue"), # gives the legend lines the correct color and width
       seg.len = 4,
       bg = "white"
)

plot(sqrt(filter.res$Sigma_tt),type="l",main=expression(X[t]),lwd=2,xlab="",ylab="",col="red",
     ylim=c(min(sqrt(smoother.res$S_smooth)),max(sqrt(filter.res$Sigma_tt))))
lines(sqrt(smoother.res$S_smooth),lwd=2,col="blue")

legend("topright", # places a legend at the appropriate place c("Health","Defense"), # puts text in the legend 
       c(
         expression(paste("Filtering standard error ",sqrt(Var(X[t*"|"*t])),sep="")),
         expression(paste("Smoothing standard error ",sqrt(Var(X[t*"|"*T])),sep=""))),
       #pt.bg = c("black","red","blue"),
       #lty=c(1,1), # gives the legend appropriate symbols (lines)       
       lwd=c(2), # line width
       lty=c(1,1,3),
       #pch=c(19,19,19),
       col=c("black","red","red"), # gives the legend lines the correct color and width
       seg.len = 4,
       bg = "white"
)

# dev.off()





# FILE = "/figures/Figure_Kalman4_new.pdf"
# pdf(file=paste(getwd(),FILE,sep=""),pointsize=10,width=7, height=5)
par(mfrow=c(2,1))
par(plt=c(.1,.95,.15,.8))

plot(X,type="l",main=expression(w[t]),lwd=2,xlab="",ylab="")
lines(filter.res$r,col="red",lwd=2)
lines(smoother.res$r,col="blue",lwd=2,lty=1)

legend("topright", # places a legend at the appropriate place c("Health","Defense"), # puts text in the legend 
       c(
         expression(w[t]),
         expression(paste("Filtered values of ",w[t],sep="")),
         expression(paste("Smoothed values of ",w[t],sep=""))),
       #pt.bg = c("black","red","blue"),
       #lty=c(1,1), # gives the legend appropriate symbols (lines)       
       lwd=c(2), # line width
       lty=c(1,1,1),
       #pch=c(19,19,19),
       col=c("black","red","blue"), # gives the legend lines the correct color and width
       seg.len = 4,
       bg = "white"
)

plot(sqrt(filter.res$Sigma_tt),type="l",
     main = "Conditional standard deviation of filtered/smoothed estimates",
     lwd=2,xlab="",ylab="",col="red",
     ylim=c(min(sqrt(smoother.res$S_smooth)),1.4*max(sqrt(filter.res$Sigma_tt))))
lines(sqrt(smoother.res$S_smooth),lwd=2,col="blue")

legend("topleft", # places a legend at the appropriate place c("Health","Defense"), # puts text in the legend 
       c(
         expression(paste("Filtering standard error ",sqrt(Var(w[t*"|"*t])),sep="")),
         expression(paste("Smoothing standard error ",sqrt(Var(w[t*"|"*T])),sep=""))),
       #pt.bg = c("black","red","blue"),
       #lty=c(1,1), # gives the legend appropriate symbols (lines)       
       lwd=c(2), # line width
       lty=c(1,1,3),
       #pch=c(19,19,19),
       col=c("red","blue"), # gives the legend lines the correct color and width
       seg.len = 4,
       bg = "white"
)

# dev.off()






