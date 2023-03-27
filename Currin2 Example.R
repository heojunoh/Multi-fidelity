install.packages("devtools")
library(devtools)
install_github("cran/MuFiCokriging")
library(MuFiCokriging)
library(lhs)
library(laGP)
source("GP.R")
source("KOH.R")
source("closed.R")
source("score.R")

fit.KOH <- function(X1, X2, X3, Y1, Y2, Y3, g=eps){ # need to change function for another example
  
  ### KOH method ###
  # Y2d3 <- fl(X3, l=3)
  # Y1d2 <- fl(X2, l=1)
  Y2d3 <- apply(X3,1,Currin2, l=2)
  Y1d2 <- apply(X2,1,Currin2, l=1)
  
  ### estimating first order ###
  fit.KOHGP1 <- KOHGP(X1, Y1)
  b1 <- 1/fit.KOHGP1$theta
  sig2_1 <- fit.KOHGP1$tau2hat
  
  ### estimating second order ###
  # KOH(X2, Y2, Y1d2)
  rho1 <- KOH(X2, Y2, Y1d2)$rho
  b2 <- 1/KOH(X2, Y2, Y1d2)$theta
  sig2_2 <- KOH(X2, Y2, Y1d2)$tau2hat
  
  ### estimating third order ###
  # KOH(X3, Y3, Y2d3)
  rho2 <- KOH(X3, Y3, Y2d3)$rho
  b3 <- 1/KOH(X3, Y3, Y2d3)$theta
  sig2_3 <- KOH(X3, Y3, Y2d3)$tau2hat
  
  return(list(b=c(b1, b2, b3), rho=c(rho1, rho2), tau2hat=c(sig2_1, sig2_2, sig2_3), g=g, X1=X1, X2=X2, X3=X3, Y1=Y1, Y2=Y2, Y3=Y3))
}

pred.KOH <- function(fit, x){ # need to change function for another example
  
  X1 <- fit$X1
  X2 <- fit$X2
  X3 <- fit$X3
  Y1 <- fit$Y1
  Y2 <- fit$Y2
  Y3 <- fit$Y3
  
  b <- fit$b
  rho <- fit$rho
  tau2hat <- fit$tau2hat
  g <- fit$g
  
  
  ### prediction of 2nd order KOH ###
  tx2 <- cbind(rho[1]*rho[2]*tau2hat[1]*covar.sep(x, X1, d=1/b[1:2], g=g), 
               rho[1]^2*rho[2]*tau2hat[1]*covar.sep(x, X2, d=1/b[1:2], g=g) + rho[2]*tau2hat[2]*covar.sep(x, X2, d=1/b[3:4], g=g),
               rho[1]^2*rho[2]^2*tau2hat[1]*covar.sep(x, X3, d=1/b[1:2], g=g) + rho[2]^2*tau2hat[2]*covar.sep(x, X3, d=1/b[3:4], g=g) + tau2hat[3]*covar.sep(x, X3, d=1/b[5:6], g=g))
  
  V1 <- tau2hat[1]*covar.sep(X1, d=1/b[1:2], g=g)
  V12 <- rho[1]*tau2hat[1]*covar.sep(X1, X2, d=1/b[1:2], g=0)
  V13 <- rho[1]*rho[2]*tau2hat[1]*covar.sep(X1, X3, d=1/b[1:2], g=0)
  V2 <- rho[1]^2*tau2hat[1]*covar.sep(X2, d=1/b[1:2], g=g) + tau2hat[2]*covar.sep(X2, d=1/b[3:4], g=g)
  V23 <- rho[1]^2*rho[2]*tau2hat[1]*covar.sep(X2, X3, d=1/b[1:2], g=0) + rho[2]*tau2hat[2]*covar.sep(X2, X3, d=1/b[3:4], g=0)
  V3 <- rho[1]^2*rho[2]^2*tau2hat[1]*covar.sep(X3, d=1/b[1:2], g=g) + rho[2]^2*tau2hat[2]*covar.sep(X3, d=1/b[3:4], g=g) + tau2hat[3]*covar.sep(X3, d=1/b[5:6], g=g)
  
  V_3 <- rbind(cbind(V1, V12, V13), cbind(t(V12), V2, V23), cbind(t(V13), t(V23), V3))
  
  mx2 <- tx2 %*% solve(V_3) %*% c(Y1, Y2, Y3)
  
  ### posterior variance ###
  koh.var2 <- pmax(0, diag(tau2hat[3]*covar.sep(as.matrix(x), d=1/b[5:6], g=g) + tau2hat[2]*rho[2]^2*covar.sep(as.matrix(x), d=1/b[3:4], g=g) + tau2hat[1]*rho[1]^2*rho[2]^2*covar.sep(as.matrix(x), d=1/b[5:6], g=g) - tx2 %*% solve(V_3+diag(g, nrow(V_3)))%*%t(tx2)))
  
  return(list(mu=mx2, sig2=koh.var2))
}

### synthetic function ###
Currin2 <- function(xx, l){
  x1 <- xx[1]
  x2 <- xx[2]
  
  term1 <- (1-exp(-1/(2*x2)))*(2300*x1^3+1900*x1^2+2092*x1+60)/(100*x1^3+500*x1^2+4*x1+20)
  term2 <- 16*2^(-l)
  
  term1 + term2*exp(-1.4*x1)*cos(3.5*pi*x2) 
}

### training data ###
n1 <- 14; n2 <- 11; n3 <- 7
set.seed(1)
X1 <- maximinLHS(n1, 2)
X2 <- maximinLHS(n2, 2)
X3 <- maximinLHS(n3, 2)

NestDesign <- NestedDesignBuild(design = list(X1,X2,X3))

X1 <- NestDesign$PX
X2 <- ExtractNestDesign(NestDesign,2)
X3 <- ExtractNestDesign(NestDesign,3)

y3 <- apply(X3,1,Currin2, l=3)
y2 <- apply(X2,1,Currin2, l=2)
y1 <- apply(X1,1,Currin2, l=1)


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


### test data ###
x <- maximinLHS(100, 2)


### closed ###
predy <- closed2(x, fit.GP1, fit.GP2new, fit.GP3new)$mu
predsig2 <- closed2(x, fit.GP1, fit.GP2new, fit.GP3new)$sig2


### KOH method ###
fit.KOH3 <- fit.KOH(X1, X2, X3, y1, y2, y3)
pred.KOH3 <- pred.KOH(fit.KOH3, x)
mx2 <- pred.KOH3$mu
koh.var2 <- pred.KOH3$sig2


### direct fitting; not using closed form. f1(u) and f_M(u) from (u, f_M(u, f1(u))) are random variables.
# w1.x <- rnorm(length(x), mean=pred.GP(fit.GP1, x)$mu, sd=sqrt(pred.GP(fit.GP1, x)$sig2)) # sample f1(x)
# w2.x <- rnorm(length(x), mean=pred.GP(fit.GP2new, cbind(x, w1.x))$mu, sd=sqrt(pred.GP(fit.GP2new, cbind(x, w1.x))$sig2))
# xxnew <- cbind(x, w2.x)
# pred3new <- pred.GP(fit.GP3new, xxnew) # not closed form
w1.x <- c(rep(NA, nrow(x)))
w2.x <- c(rep(NA, nrow(x)))
for(j in 1:nrow(x)){
  w1.x[j] <- mean(rnorm(10000, mean=pred.GP(fit.GP1, x)$mu[j], sd=sqrt(pred.GP(fit.GP1, x)$sig2[j])))
}
for(j in 1:nrow(x)){
  w2.x[j] <- mean(rnorm(10000, mean=pred.GP(fit.GP2new, cbind(x, w1.x))$mu[j], sd=sqrt(pred.GP(fit.GP2new, cbind(x, w1.x))$sig2[j])))
}

xxnew <- cbind(x, w2.x)
pred3new <- pred.GP(fit.GP3new, xxnew) # not closed form

### prediction of original GP with single fidelity ###
fit.GP3 <- GP(X3, y3)
pred3 <- pred.GP(fit.GP3, x)

### Cokriging ###
fit.muficokm <- MuFicokm(formula = list(~1,~1,~1), MuFidesign = NestDesign, #covtype="gauss",
                         response = list(y1,y2,y3), nlevel = 3)
pred.muficokm <- predict(fit.muficokm, x, "UK")


### RMSE ###
sqrt(mean((pred3$mu-apply(x,1,Currin2, l=3))^2)) # single fidelity
sqrt(mean((predy-apply(x,1,Currin2, l=3))^2)) # closed form
sqrt(mean((pred3new$mu-apply(x,1,Currin2, l=3))^2)) # not closed form
sqrt(mean((pred.muficokm$mean-apply(x,1,Currin2, l=3))^2)) # Cokriging
sqrt(mean((mx2-apply(x,1,Currin2, l=3))^2)) # KOH

mean(score(apply(x,1,Currin2, l=3), pred3$mu, pred3$sig2)) # single fidelity
mean(score(apply(x,1,Currin2, l=3), predy, predsig2)) # closed form
mean(score(apply(x,1,Currin2, l=3), pred3new$mu, pred3new$sig2)) # not closed form
mean(score(apply(x,1,Currin2, l=3), mx2, koh.var2)) # KOH

median(score(apply(x,1,Currin2, l=3), pred3$mu, pred3$sig2)) # single fidelity
median(score(apply(x,1,Currin2, l=3), predy, predsig2)) # closed form
median(score(apply(x,1,Currin2, l=3), pred3new$mu, pred3new$sig2)) # not closed form
median(score(apply(x,1,Currin2, l=3), mx2, koh.var2)) # KOH

mean(crps(apply(x,1,Currin2, l=3), pred3$mu, pred3$sig2)) # single fidelity
mean(crps(apply(x,1,Currin2, l=3), predy, predsig2)) # closed form
mean(crps(apply(x,1,Currin2, l=3), pred3new$mu, pred3new$sig2)) # not closed form
mean(crps(apply(x,1,Currin2, l=3), mx2, koh.var2)) # KOH

median(crps(apply(x,1,Currin2, l=3), pred3$mu, pred3$sig2)) # single fidelity
median(crps(apply(x,1,Currin2, l=3), predy, predsig2)) # closed form
median(crps(apply(x,1,Currin2, l=3), pred3new$mu, pred3new$sig2)) # not closed form
median(crps(apply(x,1,Currin2, l=3), mx2, koh.var2)) # KOH


sum(predsig2)
sum(pred3new$sig2)
sum(koh.var2)

