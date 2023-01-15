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
n1 <- 30; n2 <- 20; n3 <- 8 
# X1 <- maximinLHS(n1, 1) # x^L
# y1 <- fl(X1, l=1)
# X2 <- matrix(sample(X1, n2)) # x^M
# y2 <- fl(X2, l=3)
# X3 <- matrix(sample(X2, n3)) # x^H
# y3 <- fl(X3, l=5)

rep <- 100
result.3layers.rmse <- matrix(NA, rep, 3)
colnames(result.3layers.rmse) <- c("single", "closed", "direct")


for(i in 1:rep) {
  set.seed(i)

X3 <- maximinLHS(n3, 1) # x^H
y3 <- fl(X3, l=5)
X2 <- matrix(c(X3, maximinLHS(n2-n3, 1))) # x^M
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

### model fitting using (x3, f2(x3, f1(x3))) ###
# w1.x3 <- rnorm(n3, mean=pred.GP(fit.GP1, X3)$mu, sd=sqrt(pred.GP(fit.GP1, X3)$sig2)) # sample f1(x3)
# w2.x3 <- rnorm(n3, mean=pred.GP(fit.GP2new, cbind(X3, w1.x3))$mu, sd=sqrt(pred.GP(fit.GP2new, cbind(X3, w1.x3))$sig2)) # sample f2(x3, f1(x3))
w1.x3 <- pred.GP(fit.GP1, X3)$mu # can interpolate; nested
w2.x3 <- pred.GP(fit.GP2new, cbind(X3, w1.x3))$mu # can interpolate; nested
X3new <- cbind(X3, w2.x3) # combine (X3, f2(x3, f1(x3)))
fit.GP3new <- GP(X3new, y3) # model fitting for f_H(X3, f2(x3, f1(x3)))

### test data ###
x <- seq(0,1,0.001)
w1.x <- rnorm(length(x), mean=pred.GP(fit.GP1, x)$mu, sd=sqrt(pred.GP(fit.GP1, x)$sig2)) # sample f1(x)
w2.x <- rnorm(length(x), mean=pred.GP(fit.GP2new, cbind(x, w1.x))$mu, sd=sqrt(pred.GP(fit.GP2new, cbind(x, w1.x))$sig2))
xnew <- cbind(x, w2.x)

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

closed2 <- function(x, fit1, fit2, fit3){ # fit1; first layer, fit2; second layer f_M(X2, f1(X2)), fit3; last layer f_H(X3, f_M(X3, f1(X3)))
  
  d <- ncol(fit1$X)
  x <- matrix(x, ncol=d)
  
  ### second layer's output of test data 
  closed1 <- closed(x, fit1, fit2)
  x.mu <- closed1$mu
  sig2 <- closed1$sig2
  
  ### combine input ###
  # w1.x <- rnorm(length(x), mean=pred.GP(fit1, x)$mu, sd=sqrt(pred.GP(fit1, x)$sig2)) # sample f1(u)
  # x.mu <- pred.GP(fit2, cbind(x, w1.x))$mu # mean of f*(u)
  
  ### calculate the closed form ###
  X3 <- matrix(fit3$X[,-(d+1)], ncol=d)
  w2.x3 <- fit3$X[,d+1]
  y3 <- fit3$y
  n <- length(y3)
  theta <- fit3$theta
  tau2hat <- fit3$tau2hat
  
  Ci <- fit3$Ki
  a <- Ci %*% y3
  # sig2 <- pred.GP(fit2, cbind(x, w1.x))$sig2 # variance of f_M(x)
  
  # mean
  predy <- c(rep(0, length(x)))
  
  for(j in 1: length(x)){
    predv <- c(rep(1,n))
    for(i in 1:n){
      for(m in 1:d){
        predv[i] <- predv[i] * exp(-(x[j,m]-X3[i,m])^2/theta[m])
      }
      predv[i] <- predv[i] * 1/sqrt(1+2*sig2[j]/theta[d+1]) *
        exp(-(w2.x3[i]-x.mu[j])^2/(theta[d+1]+2*sig2[j]))
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
          mat[k,l] <- mat[k,l] * exp(-((x[i,m]-X3[k,m])^2+(x[i,m]-X3[l,m])^2)/theta[m]) 
        }
        mat[k,l] <- mat[k,l] * (a[k]*a[l] - tau2hat*Ci[k,l]) *
          1/sqrt(1+4*sig2[i]/theta[d+1]) * 
          exp(-((w2.x3[k]+w2.x3[l])/2-x.mu[i])^2/(theta[d+1]/2+2*sig2[i])) * 
          exp(-(w2.x3[k]-w2.x3[l])^2/(2*theta[d+1]))
      }}
    
    predsig2[i] <- pmax(0, tau2hat - predy[i]^2 + sum(mat))
  }
  
  return(list(mu=predy, sig2=predsig2))
}

predy <- closed2(x, fit.GP1, fit.GP2new, fit.GP3new)$mu
predsig2 <- closed2(x, fit.GP1, fit.GP2new, fit.GP3new)$sig2

### compared to single fidelity ###
fit.GP3 <- GP(X3, y3)
pred3 <- pred.GP(fit.GP3, x)

### direct fitting; not using closed form. f1(u) and f_M(u) from (u, f_M(u, f1(u))) are random variables.
w1.x <- rnorm(length(x), mean=pred.GP(fit.GP1, x)$mu, sd=sqrt(pred.GP(fit.GP1, x)$sig2)) # sample f1(x)
w2.x <- rnorm(length(x), mean=pred.GP(fit.GP2new, cbind(x, w1.x))$mu, sd=sqrt(pred.GP(fit.GP2new, cbind(x, w1.x))$sig2))
xnew <- cbind(x, w2.x)
pred3new <- pred.GP(fit.GP3new, xnew) # not closed form

### plot ###
# plot(x, predy, type="l", lwd=2, col=3, ylim=range(c(predy+1.96*sqrt(predsig2*length(y3)/(length(y3)-2)),predy-1.96*sqrt(predsig2*length(y3)/(length(y3)-2))))) # Green; Closed form
# lines(x, predy+1.96*sqrt(predsig2*length(y3)/(length(y3)-2)), col=3, lty=2)
# lines(x, predy-1.96*sqrt(predsig2*length(y3)/(length(y3)-2)), col=3, lty=2)
# 
# lines(x, pred3$mu, lwd=2, col=4) # Blue; Single fidelity
# lines(x, pred3$mu+1.96*sqrt(pred3$sig2*length(y3)/(length(y3)-2)), col=4, lty=2)
# lines(x, pred3$mu-1.96*sqrt(pred3$sig2*length(y3)/(length(y3)-2)), col=4, lty=2)
# 
# lines(x, pred3new$mu, lwd=2, col=7) # Yellow; Direct fitting
# lines(x, pred3new$mu+1.96*sqrt(pred3new$sig2*length(y3)/(length(y3)-2)), col=7, lty=2)
# lines(x, pred3new$mu-1.96*sqrt(pred3new$sig2*length(y3)/(length(y3)-2)), col=7, lty=2)
# 
# curve(fl(x,l=5),add=TRUE, col=1,lwd=2,lty=2) # high fidelity(TRUE); Black

### RMSE ###
result.3layers.rmse[i,1] <- sqrt(mean((pred3$mu-fl(x, l=5))^2)) # single fidelity
result.3layers.rmse[i,2] <- sqrt(mean((predy-fl(x, l=5))^2)) # closed form
result.3layers.rmse[i,3] <- sqrt(mean((pred3new$mu-fl(x, l=5))^2)) # not closed form

}

par(mfrow=c(1,1))
#RMSE comparison#
apply(result.3layers.rmse, 2, mean)
table(apply(result.3layers.rmse, 1, which.min))
boxplot(result.3layers.rmse)


