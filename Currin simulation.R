### Currin Example ###
library(lhs)
library(laGP)
source("GP.R")
source("KOH.R")
source("closed.R")
source("score.R")

### synthetic function ###
curretal88exp <- function(xx)
{
  x1 <- xx[1]
  x2 <- xx[2]
  
  fact1 <- 1 - exp(-1/(2*x2))
  fact2 <- 2300*x1^3 + 1900*x1^2 + 2092*x1 + 60
  fact3 <- 100*x1^3 + 500*x1^2 + 4*x1 + 20
  
  y <- fact1 * fact2/fact3
  return(y)
}

curretal88explc <- function(xx)
{
  x1 <- xx[1]
  x2 <- xx[2]
  
  maxarg <- max(c(0, x2-1/20))
  
  yh1 <- curretal88exp(c(x1+1/20, x2+1/20))
  yh2 <- curretal88exp(c(x1+1/20, maxarg))
  yh3 <- curretal88exp(c(x1-1/20, x2+1/20))
  yh4 <- curretal88exp(c(x1-1/20, maxarg))
  
  y <- (yh1 + yh2 + yh3 + yh4) / 4
  return(y)
}

### training data ###
n1 <- 20; n2 <- 10
d <- 2

rep <- 100
result.currin.rmse <- matrix(NA, rep, 4)
colnames(result.currin.rmse) <- c("single", "closed", "direct", "KOH")


for(i in 1:rep) {
  set.seed(i)

X2 <- maximinLHS(n2, d) # x^H
y2 <- apply(X2,1,curretal88exp)

X1 <- rbind(X2, maximinLHS(n1-n2, d)) # x^L
y1 <- apply(X1,1,curretal88explc)

### model fitting for f1 ###
eps <- sqrt(.Machine$double.eps)
fit.GP1 <- GP(X1, y1)

### model fitting using (x2, f1(x2)) ###
w1.x2 <- pred.GP(fit.GP1, X2)$mu # can interpolate; nested
X2new <- cbind(X2, w1.x2) # combine (X2, f1(x2))
fit.GP2new <- GP(X2new, y2) # model fitting for f_M(X2, f1(x2))


### test data ###
x <- maximinLHS(100, d)


### closed ###
predy <- closed(x, fit.GP1, fit.GP2new)$mu
predsig2 <- closed(x, fit.GP1, fit.GP2new)$sig2

### compared to single fidelity ###
fit.GP2 <- GP(X2, y2)
pred2 <- pred.GP(fit.GP2, x)

### direct fitting; not using closed form. f1(u) from (u, f1(u)) is random variable.
x1.mu <- rnorm(nrow(x), mean=pred.GP(fit.GP1, x)$mu, sd=sqrt(pred.GP(fit.GP1, x)$sig2)) 
xnew <- cbind(x, x1.mu) # Use mu of the input in the closed form
pred2new <- pred.GP(fit.GP2new, xnew) # not closed form

### KOH method ###
y1d2 <- apply(X2,1,curretal88explc)

### estimating first order ###
fit.KOHGP1 <- KOHGP(X1, y1)
b1 <- 1/fit.KOHGP1$theta
sig2_1 <- fit.KOHGP1$tau2hat

### estimating second order ###
# KOH(X2, y2, y1d2)
rho1 <- KOH(X2, y2, y1d2)$rho
b2 <- 1/KOH(X2, y2, y1d2)$theta
sig2_2 <- KOH(X2, y2, y1d2)$tau2hat

### prediction of 2nd order KOH ###
tx1 <- cbind(rho1*sig2_1*covar.sep(x, X1, d=1/b1, g=eps), 
             rho1^2*sig2_1*covar.sep(x, X2, d=1/b1, g=eps) + sig2_2*covar.sep(x, X2, d=1/b2, g=eps))

V1 <- sig2_1*covar.sep(X1, d=1/b1, g=eps)
V12 <- rho1*sig2_1*covar.sep(X1, X2, d=1/b1, g=0)
V2 <- rho1^2*sig2_1*covar.sep(X2, d=1/b1, g=eps) + sig2_2*covar.sep(X2, d=1/b2, g=eps)

V_2 <- rbind(cbind(V1, V12), cbind(t(V12), V2))

mx1 <- tx1 %*% solve(V_2) %*% c(y1, y2)

### posterior variance ###
koh.var1 <- pmax(0, diag(sig2_2*covar.sep(x, d=1/b2, g=eps) + sig2_1*rho1^2*covar.sep(x, d=1/b1, g=eps) - tx1 %*% solve(V_2)%*%t(tx1)))


### RMSE ###
result.currin.rmse[i,1] <- sqrt(mean((pred2$mu-apply(x,1,curretal88exp))^2)) # single fidelity
result.currin.rmse[i,2] <- sqrt(mean((predy-apply(x,1,curretal88exp))^2)) # closed form
result.currin.rmse[i,3] <- sqrt(mean((pred2new$mu-apply(x,1,curretal88exp))^2)) # not closed form
result.currin.rmse[i,4] <- sqrt(mean((mx1-apply(x,1,curretal88exp))^2)) # KOH

}

par(mfrow=c(1,1))
#RMSE comparison#
apply(result.currin.rmse, 2, mean)
table(apply(result.currin.rmse, 1, which.min))
boxplot(result.currin.rmse)


