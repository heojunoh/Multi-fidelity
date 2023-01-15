library(lhs)
library(laGP)
source("GP.R")

### synthetic function ###
fl <- function(x, l){
  term1 <- sin(2*pi*x)
  term2 <- 0.2 * sin(8*pi*x)
  
  term1 + term2*0.8^l*5 # term1 + error term
}

### training data ###
n1 <- 20; n2 <- 10;
# X1 <- maximinLHS(n1, 1) # x^L
# y1 <- fl(X1, l=1)
# X2 <- matrix(sample(X1, n2)) # x^M; nested
# y2 <- fl(X2, l=3)
X2 <- maximinLHS(n2, 1) # x^M
y2 <- fl(X2, l=3)
X1 <- matrix(c(X2, maximinLHS(n1-n2, 1))) # x^L
y1 <- fl(X1, l=1)

### model fitting for f1 ###
eps <- sqrt(.Machine$double.eps)
fit.GP1 <- GP(X1, y1)
# pred1 <- pred.GP(fit.GP1, x)

### model fitting using (x2, f1(x2)) ###
# w1.x2 <- rnorm(n2, mean=pred.GP(fit.GP1, X2)$mu, sd=sqrt(pred.GP(fit.GP1, X2)$sig2)) # sample f1(x2)
w1.x2 <- pred.GP(fit.GP1, X2)$mu # can interpolate; nested
X2new <- cbind(X2, w1.x2) # combine (X2, f1(x2))
fit.GP2new <- GP(X2new, y2) # model fitting for f_M(X2, f1(x2))

### test data ###
x <- seq(0,1,0.01)

### Function to calculate the closed form ###
closed <- function(x, fit1, fit2){ # fit1; first layer, fit2; last layer's emulator such as f_M(X2, f1(x2))
  
  d <- ncol(fit1$X)
  x <- matrix(x, ncol=d)
  x.mu <- pred.GP(fit1, x)$mu # mean of f1(u)
  
  ### calculate the closed form ###
  X2 <- matrix(fit2$X[,-(d+1)], ncol=d)
  w1.x2 <- fit2$X[,d+1]
  y2 <- fit2$y
  n <- length(y2)
  theta <- fit2$theta
  tau2hat <- fit2$tau2hat
  
  Ci <- fit2$Ki
  a <- Ci %*% y2
  sig2 <- pred.GP(fit1, x)$sig2 # variance of f1(x)
  
  # mean
  predy <- c(rep(0, length(x)))
  
  for(j in 1: length(x)){
    predv <- c(rep(1,n))
    for(i in 1:n){
      for(m in 1:d){
        predv[i] <- predv[i] * exp(-(x[j,m]-X2[i,m])^2/theta[m])
      }
      predv[i] <- predv[i] * 1/sqrt(1+2*sig2[j]/theta[d+1]) *
        exp(-(w1.x2[i]-x.mu[j])^2/(theta[d+1]+2*sig2[j]))
    }
    predy[j] <- drop(predv%*%a)
  }
  
  # var
  predsig2 <- c(rep(0, length(x)))
  
  for(i in 1: length(x)){
    mat <- matrix(1, n, n)
    for(k in 1:n){
      for(l in 1:n){
        for(m in 1:d){
          mat[k,l] <- mat[k,l] * exp(-((x[i,m]-X2[k,m])^2+(x[i,m]-X2[l,m])^2)/theta[m]) 
        }
        mat[k,l] <- mat[k,l] * (a[k]*a[l] - tau2hat*Ci[k,l]) *
          1/sqrt(1+4*sig2[i]/theta[d+1]) * 
          exp(-((w1.x2[k]+w1.x2[l])/2-x.mu[i])^2/(theta[d+1]/2+2*sig2[i])) * 
          exp(-(w1.x2[k]-w1.x2[l])^2/(2*theta[d+1]))
      }}
    
    predsig2[i] <- pmax(0, tau2hat - predy[i]^2 + sum(mat))
  }
  
  return(list(mu=predy, sig2=predsig2))
}

predy <- closed(x, fit.GP1, fit.GP2new)$mu
predsig2 <- closed(x, fit.GP1, fit.GP2new)$sig2

### compared to single fidelity ###
fit.GP2 <- GP(X2, y2)
pred2 <- pred.GP(fit.GP2, x)

### direct fitting; not using closed form. f1(u) from (u, f1(u)) is random variable.
x1.mu <- rnorm(length(x), mean=pred.GP(fit.GP1, x)$mu, sd=sqrt(pred.GP(fit.GP1, x)$sig2)) 
xnew <- cbind(x, x1.mu) # Use mu of the input in the closed form
pred2new <- pred.GP(fit.GP2new, xnew) # not closed form

### plot ###
plot(x, predy, type="l", lwd=2, col=3, ylim=range(c(predy+1.96*sqrt(predsig2*length(y2)/(length(y2)-2)),predy-1.96*sqrt(predsig2*length(y2)/(length(y2)-2))))) # Green; Closed form
lines(x, predy+1.96*sqrt(predsig2*length(y2)/(length(y2)-2)), col=3, lty=2)
lines(x, predy-1.96*sqrt(predsig2*length(y2)/(length(y2)-2)), col=3, lty=2)

lines(x, pred2$mu, lwd=2, col=4) # Blue; Single fidelity
lines(x, pred2$mu+1.96*sqrt(pred2$sig2*length(y2)/(length(y2)-2)), col=4, lty=2)
lines(x, pred2$mu-1.96*sqrt(pred2$sig2*length(y2)/(length(y2)-2)), col=4, lty=2)

lines(x, pred2new$mu, lwd=2, col=7) # Yellow; Direct fitting
lines(x, pred2new$mu+1.96*sqrt(pred2new$sig2*length(y2)/(length(y2)-2)), col=7, lty=2)
lines(x, pred2new$mu-1.96*sqrt(pred2new$sig2*length(y2)/(length(y2)-2)), col=7, lty=2)

curve(fl(x,l=3),add=TRUE, col=1,lwd=2,lty=2) # medium fidelity(TRUE); Black

### RMSE ###
sqrt(mean((predy-fl(x, l=3))^2)) # closed form
sqrt(mean((pred2$mu-fl(x, l=3))^2)) # single fidelity
sqrt(mean((pred2new$mu-fl(x, l=3))^2)) # not closed form
sum(predy-pred2new$mu) # difference between closed form and direct way

