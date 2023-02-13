library(lhs)
library(laGP)
source("GP.R")
source("KOH.R")
source("closed.R")
source("score.R")

### synthetic function ###
fl <- function(x, l){
  term1 <- sin(2*pi*x)
  term2 <- 0.2 * sin(8*pi*x)
  
  term1 + term2*0.8^l*5 + (term1+term2)^3 + exp(-2*term1*term2) #* (term1+term2)^3 # term1 + error term + interaction term
}

### training data ###
n1 <- 15; n2 <- 13; n3 <- 11

rep <- 100
result.3layers.rmse <- matrix(NA, rep, 4)
result.3layers.meanscore <- matrix(NA, rep, 4)
result.3layers.medscore <- matrix(NA, rep, 4)
result.3layers.meancrps <- matrix(NA, rep, 4)
result.3layers.medcrps <- matrix(NA, rep, 4)
colnames(result.3layers.rmse) <- c("single", "closed", "direct", "KOH")
colnames(result.3layers.meanscore) <- c("single", "closed", "direct", "KOH") # The larger, the better
colnames(result.3layers.medscore) <- c("single", "closed", "direct", "KOH") # The larger, the better
colnames(result.3layers.meancrps) <- c("single", "closed", "direct", "KOH") # The smaller, the better
colnames(result.3layers.medcrps) <- c("single", "closed", "direct", "KOH") # The smaller, the better


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
  x <- seq(0,1,0.01)
  
  
  ### closed ###
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
  koh.var1 <- pmax(0, diag(sig2_2*covar.sep(matrix(x), d=1/b2, g=eps) + sig2_1*rho1^2*covar.sep(matrix(x), d=1/b1, g=eps) - tx1 %*% solve(V_2+diag(eps, nrow(V_2)))%*%t(tx1)))
  
  ### estimating third order ###
  # KOH(X3, y3, y2d3)
  rho2 <- KOH(X3, y3, y2d3)$rho
  b3 <- 1/KOH(X3, y3, y2d3)$theta
  sig2_3 <- KOH(X3, y3, y2d3)$tau2hat
  
  ### prediction of 2nd order KOH ###
  tx2 <- cbind(rho1*rho2*sig2_1*covar.sep(x, X1, d=1/b1, g=eps), 
               rho1^2*rho2*sig2_1*covar.sep(x, X2, d=1/b1, g=eps) + rho2*sig2_2*covar.sep(x, X2, d=1/b2, g=eps),
               rho1^2*rho2^2*sig2_1*covar.sep(x, X3, d=1/b1, g=eps) + rho2^2*sig2_2*covar.sep(x, X3, d=1/b2, g=eps) + sig2_3*covar.sep(x, X3, d=1/b3, g=eps))
  
  V1 <- sig2_1*covar.sep(X1, d=1/b1, g=eps)
  V12 <- rho1*sig2_1*covar.sep(X1, X2, d=1/b1, g=0)
  V13 <- rho1*rho2*sig2_1*covar.sep(X1, X3, d=1/b1, g=0)
  V2 <- rho1^2*sig2_1*covar.sep(X2, d=1/b1, g=eps) + sig2_2*covar.sep(X2, d=1/b2, g=eps)
  V23 <- rho1^2*rho2*sig2_1*covar.sep(X2, X3, d=1/b1, g=0) + rho2*sig2_2*covar.sep(X2, X3, d=1/b2, g=0)
  V3 <- rho1^2*rho2^2*sig2_1*covar.sep(X3, d=1/b1, g=eps) + rho2^2*sig2_2*covar.sep(X3, d=1/b2, g=eps) + sig2_3*covar.sep(X3, d=1/b3, g=eps)
  
  V_3 <- rbind(cbind(V1, V12, V13), cbind(t(V12), V2, V23), cbind(t(V13), t(V23), V3))
  
  mx2 <- tx2 %*% solve(V_3) %*% c(y1, y2, y3)
  
  ### posterior variance ###
  koh.var2 <- pmax(0, diag(sig2_3*covar.sep(matrix(x), d=1/b3, g=eps) + sig2_2*rho2^2*covar.sep(matrix(x), d=1/b2, g=eps) + sig2_1*rho1^2*rho2^2*covar.sep(matrix(x), d=1/b3, g=eps) - tx2 %*% solve(V_3+diag(eps, nrow(V_3)))%*%t(tx2)))
  
  
  ### direct fitting; not using closed form. f1(u) and f_M(u) from (u, f_M(u, f1(u))) are random variables.
  # w1.x <- rnorm(length(x), mean=pred.GP(fit.GP1, x)$mu, sd=sqrt(pred.GP(fit.GP1, x)$sig2)) # sample f1(x)
  # w2.x <- rnorm(length(x), mean=pred.GP(fit.GP2new, cbind(x, w1.x))$mu, sd=sqrt(pred.GP(fit.GP2new, cbind(x, w1.x))$sig2))
  # xxnew <- cbind(x, w2.x)
  # pred3new <- pred.GP(fit.GP3new, xxnew) # not closed form
  w1.x <- c(rep(NA, length(x)))
  w2.x <- c(rep(NA, length(x)))
  for(j in 1:length(x)){
    w1.x[j] <- mean(rnorm(10000, mean=pred.GP(fit.GP1, x[j])$mu, sd=sqrt(pred.GP(fit.GP1, x[j])$sig2)))
  }
  for(j in 1:length(x)){
    w2.x[j] <- mean(rnorm(10000, mean=pred.GP(fit.GP2new, cbind(x[j], w1.x[j]))$mu, sd=sqrt(pred.GP(fit.GP2new, cbind(x[j], w1.x[j]))$sig2)))
  }
  
  xxnew <- cbind(x, w2.x)
  pred3new <- pred.GP(fit.GP3new, xxnew) # not closed form
  
  ### prediction of original GP with single fidelity ###
  fit.GP3 <- GP(X3, y3)
  pred3 <- pred.GP(fit.GP3, x)
  
  ### RMSE ###
  result.3layers.rmse[i,1] <- sqrt(mean((pred3$mu-fl(x, l=5))^2)) # single fidelity
  result.3layers.rmse[i,2] <- sqrt(mean((predy-fl(x, l=5))^2)) # closed form
  result.3layers.rmse[i,3] <- sqrt(mean((pred3new$mu-fl(x, l=5))^2)) # not closed form
  result.3layers.rmse[i,4] <- sqrt(mean((mx2-fl(x, l=5))^2)) # KOH
  
  result.3layers.meanscore[i,1] <- mean(score(fl(x, l=5), pred3$mu, pred3$sig2)) # single fidelity
  result.3layers.meanscore[i,2] <- mean(score(fl(x, l=5), predy, predsig2)) # closed form
  result.3layers.meanscore[i,3] <- mean(score(fl(x, l=5), pred3new$mu, pred3new$sig2)) # not closed form
  result.3layers.meanscore[i,4] <- mean(score(fl(x, l=5), mx2, koh.var2)) # KOH
  
  result.3layers.medscore[i,1] <- median(score(fl(x, l=5), pred3$mu, pred3$sig2)) # single fidelity
  result.3layers.medscore[i,2] <- median(score(fl(x, l=5), predy, predsig2)) # closed form
  result.3layers.medscore[i,3] <- median(score(fl(x, l=5), pred3new$mu, pred3new$sig2)) # not closed form
  result.3layers.medscore[i,4] <- median(score(fl(x, l=5), mx2, koh.var2)) # KOH
  
  result.3layers.meancrps[i,1] <- mean(crps(fl(x, l=5), pred3$mu, pred3$sig2)) # single fidelity
  result.3layers.meancrps[i,2] <- mean(crps(fl(x, l=5), predy, predsig2)) # closed form
  result.3layers.meancrps[i,3] <- mean(crps(fl(x, l=5), pred3new$mu, pred3new$sig2)) # not closed form
  result.3layers.meancrps[i,4] <- mean(crps(fl(x, l=5), mx2, koh.var2)) # KOH
  
  result.3layers.medcrps[i,1] <- median(crps(fl(x, l=5), pred3$mu, pred3$sig2)) # single fidelity
  result.3layers.medcrps[i,2] <- median(crps(fl(x, l=5), predy, predsig2)) # closed form
  result.3layers.medcrps[i,3] <- median(crps(fl(x, l=5), pred3new$mu, pred3new$sig2)) # not closed form
  result.3layers.medcrps[i,4] <- median(crps(fl(x, l=5), mx2, koh.var2)) # KOH
  
}

par(mfrow=c(1,1))
#RMSE comparison#
apply(result.3layers.rmse, 2, mean)
table(apply(result.3layers.rmse, 1, which.min))
boxplot(result.3layers.rmse)

#score comparison, The larger, the better
apply(result.3layers.meanscore, 2, mean)
table(apply(result.3layers.meanscore, 1, which.max))
boxplot(result.3layers.meanscore)

#score comparison, The larger, the better
apply(result.3layers.medscore, 2, mean)
table(apply(result.3layers.medscore, 1, which.max))
boxplot(result.3layers.medscore)

#CRPS comparison, The smaller, the better
apply(result.3layers.meancrps, 2, mean)
table(apply(result.3layers.meancrps, 1, which.min))
boxplot(result.3layers.meancrps)

#CRPS comparison, The smaller, the better
apply(result.3layers.medcrps, 2, mean)
table(apply(result.3layers.medcrps, 1, which.min))
boxplot(result.3layers.medcrps)

