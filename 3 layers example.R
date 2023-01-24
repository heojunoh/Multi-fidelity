library(lhs)
library(laGP)
source("GP.R")
source("KOH.R")

### synthetic function ###
fl <- function(x, l){
  term1 <- sin(2*pi*x)
  term2 <- 0.2 * sin(8*pi*x)

  term1 + term2*0.8^l*5  # term1 + error term
}
set.seed(1)
### training data ###
n1 <- 30; n2 <- 20; n3 <- 8

X3 <- maximinLHS(n3, 1) # x^H
y3 <- fl(X3, l=5)
X2 <- matrix(c(X3, maximinLHS(n2-n3, 1))) # x^M
y2 <- fl(X2, l=3)
X1 <- matrix(c(X2, maximinLHS(n1-n2, 1))) # x^L
y1 <- fl(X1, l=1)

### model fitting for f1 ###
eps <- sqrt(.Machine$double.eps)
fit.GP1 <- GP(X1, y1)

### model fitting using (x2, f1(x2)) ###
w1.x2 <- pred.GP(fit.GP1, X2)$mu # can interpolate; nested
X2new <- cbind(X2, w1.x2) # combine (X2, f1(x2))
fit.GP2new <- GP(X2new, y2) # model fitting for f_M(X2, f1(x2))

### model fitting using (x3, f2(x3, f1(x3))) ###
w1.x3 <- pred.GP(fit.GP1, X3)$mu # can interpolate; nested
w2.x3 <- pred.GP(fit.GP2new, cbind(X3, w1.x3))$mu # can interpolate; nested
X3new <- cbind(X3, w2.x3) # combine (X3, f2(x3, f1(x3)))
fit.GP3new <- GP(X3new, y3) # model fitting for f_H(X3, f2(x3, f1(x3)))

### Function to calculate the closed form ###
closed <- function(x, fit1, fit2){ # fit1; first layer, fit2; last layer's emulator such as f_M(X2, f1(x2))
  
  d <- ncol(fit1$X)
  x <- matrix(x, ncol=d)
  x.mu <- pred.GP(fit1, x)$mu # mean of f1(u)
  sig2 <- pred.GP(fit1, x)$sig2 # variance of f1(x)
  
  ### calculate the closed form ###
  X2 <- matrix(fit2$X[,-(d+1)], ncol=d)
  w1.x2 <- fit2$X[,d+1]
  y2 <- fit2$y
  n <- length(y2)
  theta <- fit2$theta
  tau2hat <- fit2$tau2hat
  
  Ci <- fit2$Ki
  a <- Ci %*% (y2 * attr(y2, "scaled:scale") + attr(y2, "scaled:center"))
  
  ### scale new inputs ###
  x <- t((t(x)-attr(fit2$X,"scaled:center")[1:d])/attr(fit2$X,"scaled:scale")[1:d])
  x.mu <- t((t(x.mu)-attr(fit2$X,"scaled:center")[d+1])/attr(fit2$X,"scaled:scale")[d+1])
  
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
  
  ### calculate the closed form ###
  X3 <- matrix(fit3$X[,-(d+1)], ncol=d)
  w2.x3 <- fit3$X[,d+1]
  y3 <- fit3$y
  n <- length(y3)
  theta <- fit3$theta
  tau2hat <- fit3$tau2hat
  
  Ci <- fit3$Ki
  a <- Ci %*% (y3 * attr(y3, "scaled:scale") + attr(y3, "scaled:center"))
  
  ### scale new inputs ###
  x <- t((t(x)-attr(fit3$X,"scaled:center")[1:d])/attr(fit3$X,"scaled:scale")[1:d])
  x.mu <- t((t(x.mu)-attr(fit3$X,"scaled:center")[d+1])/attr(fit3$X,"scaled:scale")[d+1])

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

### test data ###
x <- seq(0,1,0.001)

predy <- closed2(x, fit.GP1, fit.GP2new, fit.GP3new)$mu
predsig2 <- closed2(x, fit.GP1, fit.GP2new, fit.GP3new)$sig2


### KOH method ###
y2d3 <- fl(X3, l=3)
y1d2 <- fl(X2, l=1)

### estimating first order ###
fit.KOHGP1 <- KOHGP(X1, y1)
b1 <- 1/fit.KOHGP1$theta
sig2_1 <- fit.KOHGP1$tau2hat

### estimating second order ###
KOH(X2, y2, y1d2)
rho1 <- KOH(X2, y2, y1d2)$rho
b2 <- 1/KOH(X2, y2, y1d2)$theta
sig2_2 <- KOH(X2, y2, y1d2)$tau2hat

### prediction of 2nd order KOH ###
tx1 <- cbind(rho1*sig2_1*covar.sep(x, X1, d=1/b1, g=eps), 
             rho1^2*sig2_1*covar.sep(x, X2, d=1/b1, g=eps) + sig2_2*covar.sep(x, X2, d=1/b2, g=eps))

V1 <- sig2_1*covar.sep(X1, d=1/b1, g=eps)
V12 <- rho1*sig2_1*covar.sep(X1, X2, d=1/b1, g=eps)
V2 <- rho1^2*sig2_1*covar.sep(X2, d=1/b1, g=eps) + sig2_2*covar.sep(X2, d=1/b2, g=eps)

V_2 <- rbind(cbind(V1, V12), cbind(t(V12), V2))

mx1 <- tx1 %*% solve(V_2) %*% c(y1, y2)

### posterior variance ###
koh.var1 <- pmax(0, diag(sig2_2*covar.sep(matrix(x), d=1/b2, g=eps) + sig2_1*rho1^2*covar.sep(matrix(x), d=1/b1, g=eps) - tx1 %*% solve(V_2)%*%t(tx1)))

### estimating third order ###
KOH(X3, y3, y2d3)
rho2 <- KOH(X3, y3, y2d3)$rho
b3 <- 1/KOH(X3, y3, y2d3)$theta
sig2_3 <- KOH(X3, y3, y2d3)$tau2hat

### prediction of 2nd order KOH ###
tx2 <- cbind(rho1*rho2*sig2_1*covar.sep(x, X1, d=1/b1, g=eps), 
             rho1^2*rho2*sig2_1*covar.sep(x, X2, d=1/b1, g=eps) + rho2*sig2_2*covar.sep(x, X2, d=1/b2, g=eps),
             rho1^2*rho2^2*sig2_1*covar.sep(x, X3, d=1/b1, g=eps) + rho2^2*sig2_2*covar.sep(x, X3, d=1/b2, g=eps) + sig2_3*covar.sep(x, X3, d=1/b3, g=eps))

V1 <- sig2_1*covar.sep(X1, d=1/b1, g=eps)
V12 <- rho1*sig2_1*covar.sep(X1, X2, d=1/b1, g=eps)
V13 <- rho1*rho2*sig2_1*covar.sep(X1, X3, d=1/b1, g=eps)
V2 <- rho1^2*sig2_1*covar.sep(X2, d=1/b1, g=eps) + sig2_2*covar.sep(X2, d=1/b2, g=eps)
V23 <- rho1^2*rho2*sig2_1*covar.sep(X2, X3, d=1/b1, g=eps) + rho2*sig2_2*covar.sep(X2, X3, d=1/b2, g=eps)
V3 <- rho1^2*rho2^2*sig2_1*covar.sep(X3, d=1/b1, g=eps) + rho2^2*sig2_2*covar.sep(X3, d=1/b2, g=eps) + sig2_3*covar.sep(X3, d=1/b3, g=eps)

V_3 <- rbind(cbind(V1, V12, V13), cbind(t(V12), V2, V23), cbind(t(V13), t(V23), V3))

mx2 <- tx2 %*% solve(V_3) %*% c(y1, y2, y3)

### posterior variance ###
koh.var2 <- pmax(0, diag(sig2_3*covar.sep(matrix(x), d=1/b3, g=eps) + sig2_2*rho2^2*covar.sep(matrix(x), d=1/b2, g=eps) + sig2_1*rho1^2*rho2^2*covar.sep(matrix(x), d=1/b3, g=eps) - tx2 %*% solve(V_3)%*%t(tx2)))


### direct fitting; not using closed form. f1(u) and f_M(u) from (u, f_M(u, f1(u))) are random variables.
w1.x <- rnorm(length(x), mean=pred.GP(fit.GP1, x)$mu, sd=sqrt(pred.GP(fit.GP1, x)$sig2)) # sample f1(x)
w2.x <- rnorm(length(x), mean=pred.GP(fit.GP2new, cbind(x, w1.x))$mu, sd=sqrt(pred.GP(fit.GP2new, cbind(x, w1.x))$sig2))
xxnew <- cbind(x, w2.x)
pred3new <- pred.GP(fit.GP3new, xxnew) # not closed form

### prediction of original GP with single fidelity ###
fit.GP3 <- GP(X3, y3)
pred3 <- pred.GP(fit.GP3, x)

### plot ###
plot(x, predy, type="l", lwd=2, col=3, 
     ylim=range(c(predy+1.96*sqrt(predsig2*length(y3)/(length(y3)-2)), predy-1.96*sqrt(predsig2*length(y3)/(length(y3)-2)), 
                  pred3$mu+1.96*sqrt(pred3$sig2*length(y3)/(length(y3)-2)), pred3$mu-1.96*sqrt(pred3$sig2*length(y3)/(length(y3)-2)),
                  mx2+1.96*sqrt(koh.var2*length(y3)/(length(y3)-2)), mx2-1.96*sqrt(koh.var2*length(y3)/(length(y3)-2)),
                  pred3$mu+1.96*sqrt(pred3$sig2*length(y3)/(length(y3)-2)), pred3$mu-1.96*sqrt(pred3$sig2*length(y3)/(length(y3)-2))
     ))) # Green; Closed form
lines(x, predy+1.96*sqrt(predsig2*length(y3)/(length(y3)-2)), col=3, lty=2)
lines(x, predy-1.96*sqrt(predsig2*length(y3)/(length(y3)-2)), col=3, lty=2)

lines(x, pred3new$mu, lwd=2, col=6) # Purple; Direct fitting
lines(x, pred3new$mu+1.96*sqrt(pred3new$sig2*length(y3)/(length(y3)-2)), col=6, lty=2)
lines(x, pred3new$mu-1.96*sqrt(pred3new$sig2*length(y3)/(length(y3)-2)), col=6, lty=2)

lines(x, pred3$mu, lwd=2, col=4) # Blue; Single fidelity
lines(x, pred3$mu+1.96*sqrt(pred3$sig2*length(y3)/(length(y3)-2)), col=4, lty=2)
lines(x, pred3$mu-1.96*sqrt(pred3$sig2*length(y3)/(length(y3)-2)), col=4, lty=2)

lines(x, mx2, lwd=2, col=7) # Yellow; KOH
lines(x, mx2+1.96*sqrt(koh.var2*length(y3)/(length(y3)-2)), col=7, lty=2)
lines(x, mx2-1.96*sqrt(koh.var2*length(y3)/(length(y3)-2)), col=7, lty=2)

curve(fl(x,l=5),add=TRUE, col=1,lwd=2,lty=2) # high fidelity(TRUE); Black

### RMSE ###
sqrt(mean((predy-fl(x, l=5))^2)) # closed form
sqrt(mean((pred3new$mu-fl(x, l=5))^2)) # not closed form
sqrt(mean((pred3$mu-fl(x, l=5))^2)) # single fidelity
sqrt(mean((mx2-fl(x, l=5))^2)) # KOH


sum(predsig2)
sum(pred3new$sig2)
sum(koh.var1)

