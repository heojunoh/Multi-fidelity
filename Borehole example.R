### Borehole Example ###
library(lhs)
library(laGP)
source("GP.R")
source("KOH.R")
source("closed.R")

### synthetic function ###
borehole <- function(xx)
{
  rw <- xx[1]
  r  <- xx[2]
  Tu <- xx[3]
  Hu <- xx[4]
  Tl <- xx[5]
  Hl <- xx[6]
  L  <- xx[7]
  Kw <- xx[8]
  
  frac1 <- 2 * pi * Tu * (Hu-Hl)
  
  frac2a <- 2*L*Tu / (log(r/rw)*rw^2*Kw)
  frac2b <- Tu / Tl
  frac2 <- log(r/rw) * (1+frac2a+frac2b)
  
  y <- frac1 / frac2
  return(y)
}

boreholelow <- function(xx)
{ 
  rw <- xx[1]
  r  <- xx[2]
  Tu <- xx[3]
  Hu <- xx[4]
  Tl <- xx[5]
  Hl <- xx[6]
  L  <- xx[7]
  Kw <- xx[8]
  
  frac1 <- 5 * Tu * (Hu-Hl) #+ (Tu * Kw) # Tu * Kw is added
  
  frac2a <- 2*L*Tu / (log(r/rw)*rw^2*Kw)
  frac2b <- Tu / Tl
  frac2 <- log(r/rw) * (1.5+frac2a+frac2b)
  
  y <- frac1 / frac2
  return(y)
}
set.seed(1)
output.f <- function(x){
  factor_range <- list("rw" = c(0.05, 0.15), "r" = c(100, 50000),
                       "Tu" = c(63070, 115600), "Hu" = c(990, 1110),
                       "Tl" = c(63.1, 116), "Hl" = c(700, 820),
                       "L" = c(1120, 1680), "Kw" = c(9855, 12045))
  for(i in 1:length(factor_range)) x[i] <- factor_range[[i]][1] + x[i] * diff(factor_range[[i]])
  borehole(x[1:8])
} 

outputlow.f <- function(x){
  factor_range <- list("rw" = c(0.05, 0.15), "r" = c(100, 50000),
                       "Tu" = c(63070, 115600), "Hu" = c(990, 1110),
                       "Tl" = c(63.1, 116), "Hl" = c(700, 820),
                       "L" = c(1120, 1680), "Kw" = c(9855, 12045))
  for(i in 1:length(factor_range)) x[i] <- factor_range[[i]][1] + x[i] * diff(factor_range[[i]])
  boreholelow(x[1:8])
} 

### training data ###
n1 <- 80; n2 <- 40
d <- 8

X2 <- maximinLHS(n2, d) # x^H
y2 <- apply(X2,1,output.f)

X1 <- rbind(X2, maximinLHS(n1-n2, d)) # x^L
y1 <- apply(X1,1,outputlow.f)

### model fitting for f1 ###
eps <- sqrt(.Machine$double.eps)
fit.GP1 <- GP(X1, y1)
fit.GP1$theta

### model fitting using (x2, f1(x2)) ###
w1.x2 <- pred.GP(fit.GP1, X2)$mu # can interpolate; nested
X2new <- cbind(X2, w1.x2) # combine (X2, f1(x2)) 
fit.GP2new <- GP(X2new, y2) # model fitting for f_M(X2, f1(x2))
fit.GP2new$theta

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
y1d2 <- apply(X2,1,outputlow.f)

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

V_2 <- rbind(cbind(V1, V12), cbind(t(V12), V2))+diag(eps,nrow(X1)+nrow(X2))

mx1 <- tx1 %*% solve(V_2) %*% c(y1, y2)

### posterior variance ###
koh.var1 <- pmax(0, diag(sig2_2*covar.sep(x, d=1/b2, g=eps) + sig2_1*rho1^2*covar.sep(x, d=1/b1, g=eps) - tx1 %*% solve(V_2)%*%t(tx1)))


### RMSE ###
sqrt(mean((predy-apply(x,1,output.f))^2)) # closed form
sqrt(mean((pred2new$mu-apply(x,1,output.f))^2)) # not closed form
sqrt(mean((pred2$mu-apply(x,1,output.f))^2)) # single fidelity
sqrt(mean((mx1-apply(x,1,output.f))^2)) # KOH

mean(score(apply(x,1,output.f), pred2$mu, pred2$sig2)) # single fidelity
mean(score(apply(x,1,output.f), predy, predsig2)) # closed form
mean(score(apply(x,1,output.f), pred2new$mu, pred2new$sig2)) # not closed form
mean(score(apply(x,1,output.f), mx1, koh.var1)) # KOH

median(score(apply(x,1,output.f), pred2$mu, pred2$sig2)) # single fidelity
median(score(apply(x,1,output.f), predy, predsig2)) # closed form
median(score(apply(x,1,output.f), pred2new$mu, pred2new$sig2)) # not closed form
median(score(apply(x,1,output.f), mx1, koh.var1)) # KOH

mean(crps(apply(x,1,output.f), pred2$mu, pred2$sig2)) # single fidelity
mean(crps(apply(x,1,output.f), predy, predsig2)) # closed form
mean(crps(apply(x,1,output.f), pred2new$mu, pred2new$sig2)) # not closed form
mean(crps(apply(x,1,output.f), mx1, koh.var1)) # KOH

median(crps(apply(x,1,output.f), pred2$mu, pred2$sig2)) # single fidelity
median(crps(apply(x,1,output.f), predy, predsig2)) # closed form
median(crps(apply(x,1,output.f), pred2new$mu, pred2new$sig2)) # not closed form
median(crps(apply(x,1,output.f), mx1, koh.var1)) # KOH


sum(predsig2)
sum(pred2new$sig2)
sum(koh.var1)

