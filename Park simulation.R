### Park Example ###
library(lhs)
library(laGP)
source("GP.R")

### synthetic function ###
park91b <- function(xx)
{
  x1 <- xx[1]
  x2 <- xx[2]
  x3 <- xx[3]
  x4 <- xx[4]
  
  term1 <- (2/3) * exp(x1+x2)
  term2 <- -x4 * sin(x3)
  term3 <- x3
  
  y <- term1 + term2 + term3
  return(y)
}

park91blc <- function(xx)
{
  yh <- park91b(xx)
  
  y <- 1.2*yh -1
  return(y)
}

### training data ###
n1 <- 20; n2 <- 10
d <- 4

rep <- 100
result.park.rmse <- matrix(NA, rep, 3)
colnames(result.park.rmse) <- c("single", "closed", "direct")


for(i in 1:rep) {
  set.seed(i)

X2 <- maximinLHS(n2, d) # x^H
y2 <- apply(X2,1,park91b)

X1 <- rbind(X2, maximinLHS(n1-n2, d)) # x^L
y1 <- apply(X1,1,park91blc)

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
x <- maximinLHS(100, d)

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
  predy <- c(rep(0, nrow(x)))
  
  for(j in 1: nrow(x)){
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
  predsig2 <- c(rep(0, nrow(x)))
  
  for(i in 1: nrow(x)){
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
x1.mu <- rnorm(nrow(x), mean=pred.GP(fit.GP1, x)$mu, sd=sqrt(pred.GP(fit.GP1, x)$sig2)) 
xnew <- cbind(x, x1.mu) # Use mu of the input in the closed form
pred2new <- pred.GP(fit.GP2new, xnew) # not closed form


### RMSE ###
result.park.rmse[i,1] <- sqrt(mean((pred2$mu-apply(x,1,park91b))^2)) # single fidelity
result.park.rmse[i,2] <- sqrt(mean((predy-apply(x,1,park91b))^2)) # closed form
result.park.rmse[i,3] <- sqrt(mean((pred2new$mu-apply(x,1,park91b))^2)) # not closed form
}

par(mfrow=c(1,1))
#RMSE comparison#
apply(result.park.rmse, 2, mean)
table(apply(result.park.rmse, 1, which.min))
boxplot(result.park.rmse)



