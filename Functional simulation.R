### functional Example ###
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
  Y2d3 <- f2(X3)
  Y1d2 <- f1(X2)
  
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
  tx2 <- cbind(rho[1]*rho[2]*tau2hat[1]*covar.sep(x, X1, d=1/b[1], g=g), 
               rho[1]^2*rho[2]*tau2hat[1]*covar.sep(x, X2, d=1/b[1], g=g) + rho[2]*tau2hat[2]*covar.sep(x, X2, d=1/b[2], g=g),
               rho[1]^2*rho[2]^2*tau2hat[1]*covar.sep(x, X3, d=1/b[1], g=g) + rho[2]^2*tau2hat[2]*covar.sep(x, X3, d=1/b[2], g=g) + tau2hat[3]*covar.sep(x, X3, d=1/b[3], g=g))
  
  V1 <- tau2hat[1]*covar.sep(X1, d=1/b[1], g=g)
  V12 <- rho[1]*tau2hat[1]*covar.sep(X1, X2, d=1/b[1], g=0)
  V13 <- rho[1]*rho[2]*tau2hat[1]*covar.sep(X1, X3, d=1/b[1], g=0)
  V2 <- rho[1]^2*tau2hat[1]*covar.sep(X2, d=1/b[1], g=g) + tau2hat[2]*covar.sep(X2, d=1/b[2], g=g)
  V23 <- rho[1]^2*rho[2]*tau2hat[1]*covar.sep(X2, X3, d=1/b[1], g=0) + rho[2]*tau2hat[2]*covar.sep(X2, X3, d=1/b[2], g=0)
  V3 <- rho[1]^2*rho[2]^2*tau2hat[1]*covar.sep(X3, d=1/b[1], g=g) + rho[2]^2*tau2hat[2]*covar.sep(X3, d=1/b[2], g=g) + tau2hat[3]*covar.sep(X3, d=1/b[3], g=g)
  
  V_3 <- rbind(cbind(V1, V12, V13), cbind(t(V12), V2, V23), cbind(t(V13), t(V23), V3))
  
  mx2 <- tx2 %*% solve(V_3) %*% c(Y1, Y2, Y3)
  
  ### posterior variance ###
  koh.var2 <- pmax(0, diag(tau2hat[3]*covar.sep(as.matrix(x), d=1/b[3], g=g) + tau2hat[2]*rho[2]^2*covar.sep(as.matrix(x), d=1/b[2], g=g) + tau2hat[1]*rho[1]^2*rho[2]^2*covar.sep(as.matrix(x), d=1/b[3], g=g) - tx2 %*% solve(V_3+diag(g, nrow(V_3)))%*%t(tx2)))
  
  return(list(mu=mx2, sig2=koh.var2))
}

### synthetic function ###
fl <- function(x, l){
  term1 <- sin(2*pi*x)
  term2 <- 0.2 * sin(8*pi*x)
  
  term1 + term2*5*0.8^l + (term1+term2)^3 + exp(-2*term1*term2)
}

f1 <- function(x)
{
  fl(x, Inf)
}

f2 <- function(x)
{ 
  y1 <- f1(x)
  1/10*(exp(y1) + (cos(y1) + y1)^3)
  # y1^2 - sin(abs(y1))
  # exp(-1/2*y1) + log(1/3+abs(y1))
  # 1/8*(exp(y1) + (cos(y1) + y1)^3 + y1^2 - sin(abs(y1)) + exp(-1/2*y1) + log(1/3+abs(y1)))
}

f3 <- function(x)
{ 
  y2 <- f2(x)
  # -(sin(log(1+y2)) - cos(y2))^2 + 1/2*exp(1+sin(x*pi)) 
  - cos(y2)^2 + 1/2*exp(1+sin(x*pi))
}

x <- seq(0,1,0.01)
f1(x)
f2(x)
f3(x)

plot(x, f1(x), type="l", lwd=2, col="red", ylim=c(min(c(f1(x), f2(x), f3(x))), max(c(f1(x), f2(x), f3(x)))))
lines(x, f2(x), col="orange",lwd=2,lty=1)
curve(f3(x),add=TRUE, col="green",lwd=2,lty=1) # high fidelity(TRUE); Black

### training data ###
n1 <- 24; n2 <- 22; n3 <- 20

rep <- 100
result.functional.rmse <- matrix(NA, rep, 5)
result.functional.meanscore <- matrix(NA, rep, 5)
result.functional.medscore <- matrix(NA, rep, 5)
result.functional.meancrps <- matrix(NA, rep, 5)
result.functional.medcrps <- matrix(NA, rep, 5)
colnames(result.functional.rmse) <- c("single", "closed", "direct", "Cokriging", "KOH")
colnames(result.functional.meanscore) <- c("single", "closed", "direct", "Cokriging", "KOH") # The larger, the better
colnames(result.functional.medscore) <- c("single", "closed", "direct", "Cokriging", "KOH") # The larger, the better
colnames(result.functional.meancrps) <- c("single", "closed", "direct", "Cokriging", "KOH") # The smaller, the better
colnames(result.functional.medcrps) <- c("single", "closed", "direct", "Cokriging", "KOH") # The smaller, the better

for(i in 1:rep) {
  set.seed(i)
  
  X1 <- maximinLHS(n1, 1)
  X2 <- maximinLHS(n2, 1)
  X3 <- maximinLHS(n3, 1)
  
  NestDesign <- NestedDesignBuild(design = list(X1,X2,X3))
  
  X1 <- NestDesign$PX
  X2 <- ExtractNestDesign(NestDesign,2)
  X3 <- ExtractNestDesign(NestDesign,3)
  
  y1 <- f1(X1)
  y2 <- f2(X2)
  y3 <- f3(X3)
  
  
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
  fit.KOH3 <- fit.KOH(X1, X2, X3, y1, y2, y3)
  pred.KOH3 <- pred.KOH(fit.KOH3, x)
  mx2 <- pred.KOH3$mu
  koh.var2 <- pred.KOH3$sig2
  
  
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
  
  ### Cokriging ###
  fit.muficokm <- MuFicokm(formula = list(~1,~1,~1), MuFidesign = NestDesign, #covtype="gauss",
                           coef.trend = list(0,c(0,0),c(0,0)), response = list(y1,y2,y3), nlevel = 3)
  pred.muficokm <- predict(fit.muficokm, x, "SK")
  
  ### RMSE ###
  result.functional.rmse[i,1] <- sqrt(mean((pred3$mu-f3(x))^2)) # single fidelity
  result.functional.rmse[i,2] <- sqrt(mean((predy-f3(x))^2)) # closed form
  result.functional.rmse[i,3] <- sqrt(mean((pred3new$mu-f3(x))^2)) # not closed form
  result.functional.rmse[i,4] <- sqrt(mean((pred.muficokm$mean-f3(x))^2)) # Cokriging
  result.functional.rmse[i,5] <- sqrt(mean((mx2-f3(x))^2)) # KOH
  
  result.functional.meanscore[i,1] <- mean(score(f3(x), pred3$mu, pred3$sig2)) # single fidelity
  result.functional.meanscore[i,2] <- mean(score(f3(x), predy, predsig2)) # closed form
  result.functional.meanscore[i,3] <- mean(score(f3(x), pred3new$mu, pred3new$sig2)) # not closed form
  result.functional.meanscore[i,4] <- mean(score(f3(x), pred.muficokm$mean, pred.muficokm$sig2)) # Cokriging
  result.functional.meanscore[i,5] <- mean(score(f3(x), mx2, koh.var2)) # KOH
  
  result.functional.medscore[i,1] <- median(score(f3(x), pred3$mu, pred3$sig2)) # single fidelity
  result.functional.medscore[i,2] <- median(score(f3(x), predy, predsig2)) # closed form
  result.functional.medscore[i,3] <- median(score(f3(x), pred3new$mu, pred3new$sig2)) # not closed form
  result.functional.medscore[i,4] <- median(score(f3(x), pred.muficokm$mean, pred.muficokm$sig2)) # Cokriging
  result.functional.medscore[i,5] <- median(score(f3(x), mx2, koh.var2)) # KOH
  
  result.functional.meancrps[i,1] <- mean(crps(f3(x), pred3$mu, pred3$sig2)) # single fidelity
  result.functional.meancrps[i,2] <- mean(crps(f3(x), predy, predsig2)) # closed form
  result.functional.meancrps[i,3] <- mean(crps(f3(x), pred3new$mu, pred3new$sig2)) # not closed form
  result.functional.meancrps[i,4] <- mean(crps(f3(x), pred.muficokm$mean, pred.muficokm$sig2)) # Cokriging
  result.functional.meancrps[i,5] <- mean(crps(f3(x), mx2, koh.var2)) # KOH
  
  result.functional.medcrps[i,1] <- median(crps(f3(x), pred3$mu, pred3$sig2)) # single fidelity
  result.functional.medcrps[i,2] <- median(crps(f3(x), predy, predsig2)) # closed form
  result.functional.medcrps[i,3] <- median(crps(f3(x), pred3new$mu, pred3new$sig2)) # not closed form
  result.functional.medcrps[i,4] <- median(crps(f3(x), pred.muficokm$mean, pred.muficokm$sig2)) # Cokriging
  result.functional.medcrps[i,5] <- median(crps(f3(x), mx2, koh.var2)) # KOH
  
}


par(mfrow=c(1,1))
#RMSE comparison#
apply(result.functional.rmse, 2, mean)
table(apply(result.functional.rmse, 1, which.min))
boxplot(result.functional.rmse)

#score comparison, The larger, the better
apply(result.functional.meanscore, 2, mean)
table(apply(result.functional.meanscore, 1, which.max))
boxplot(result.functional.meanscore)

#score comparison, The larger, the better
apply(result.functional.medscore, 2, mean)
table(apply(result.functional.medscore, 1, which.max))
boxplot(result.functional.medscore)

#CRPS comparison, The smaller, the better
apply(result.functional.meancrps, 2, mean)
table(apply(result.functional.meancrps, 1, which.min))
boxplot(result.functional.meancrps)

#CRPS comparison, The smaller, the better
apply(result.functional.medcrps, 2, mean)
table(apply(result.functional.medcrps, 1, which.min))
boxplot(result.functional.medcrps)


# closedrmse <- matrix(c(0,0,0), nrow=1)
# closedmeanscore <- matrix(c(0,0,0), nrow=1)
# closedmedscore <- matrix(c(0,0,0), nrow=1)
# closedmeancrps <- matrix(c(0,0,0), nrow=1)
# closedmedcrps <- matrix(c(0,0,0), nrow=1)
# kohrmse <- matrix(c(0,0,0), nrow=1)
# kohmeanscore <- matrix(c(0,0,0), nrow=1)
# kohmedscore <- matrix(c(0,0,0), nrow=1)
# kohmeancrps <- matrix(c(0,0,0), nrow=1)
# kohmedcrps <- matrix(c(0,0,0), nrow=1)


closedrmse <- rbind(closedrmse, c(sort(result.functional.rmse[,2])[5], apply(result.functional.rmse, 2, mean)[2], sort(result.functional.rmse[,2])[96]))
closedmeanscore <- rbind(closedmeanscore, c(sort(result.functional.meanscore[,2])[5], apply(result.functional.meanscore, 2, mean)[2], sort(result.functional.meanscore[,2])[96]))
closedmedscore <- rbind(closedmedscore, c(sort(result.functional.medscore[,2])[5], apply(result.functional.medscore, 2, mean)[2], sort(result.functional.medscore[,2])[96]))
closedmeancrps <- rbind(closedmeancrps, c(sort(result.functional.meancrps[,2])[5], apply(result.functional.meancrps, 2, mean)[2], sort(result.functional.meancrps[,2])[96]))
closedmedcrps <- rbind(closedmedcrps, c(sort(result.functional.medcrps[,2])[5], apply(result.functional.medcrps, 2, mean)[2], sort(result.functional.medcrps[,2])[96]))

kohrmse <- rbind(kohrmse, c(sort(result.functional.rmse[,4])[5], apply(result.functional.rmse, 2, mean)[4], sort(result.functional.rmse[,4])[96]))
kohmeanscore <- rbind(kohmeanscore, c(sort(result.functional.meanscore[,4])[5], apply(result.functional.meanscore, 2, mean)[4], sort(result.functional.meanscore[,4])[96]))
kohmedscore <- rbind(kohmedscore, c(sort(result.functional.medscore[,4])[5], apply(result.functional.medscore, 2, mean)[4], sort(result.functional.medscore[,4])[96]))
kohmeancrps <- rbind(kohmeancrps, c(sort(result.functional.meancrps[,4])[5], apply(result.functional.meancrps, 2, mean)[4], sort(result.functional.meancrps[,4])[96]))
kohmedcrps <- rbind(kohmedcrps, c(sort(result.functional.medcrps[,4])[5], apply(result.functional.medcrps, 2, mean)[4], sort(result.functional.medcrps[,4])[96]))


# closedrmse <- closedrmse[-1,]
# closedmeanscore <- closedmeanscore[-1,]
# closedmedscore <- closedmedscore[-1,]
# closedmeancrps <- closedmeancrps[-1,]
# closedmedcrps <- closedmedcrps[-1,]
# 
# kohrmse <- kohrmse[-1,]
# kohmeanscore <- kohmeanscore[-1,]
# kohmedscore <- kohmedscore[-1,]
# kohmeancrps <- kohmeancrps[-1,]
# kohmedcrps <- kohmedcrps[-1,]


### RMSE ###
plot(c(30,40,50,60,70), closedrmse[,2], type="l", lwd=2, col=3,  # Green; closed
     ylim=c(min(c(closedrmse, kohrmse)), max(c(closedrmse, kohrmse)))) 
lines(c(30,40,50,60,70), closedrmse[,1], col=3, lty=2)
lines(c(30,40,50,60,70), closedrmse[,3], col=3, lty=2)

lines(c(30,40,50,60,70), kohrmse[,2], lwd=2, col=7) # Yellow; KOH
lines(c(30,40,50,60,70), kohrmse[,1], col=7, lty=2)
lines(c(30,40,50,60,70), kohrmse[,3], col=7, lty=2)

### Mean score ###
plot(c(30,40,50,60,70), closedmeanscore[,2], type="l", lwd=2, col=3,  # Green; closed
     ylim=c(min(c(closedmeanscore, kohmeanscore)), max(c(closedmeanscore, kohmeanscore)))) 
lines(c(30,40,50,60,70), closedmeanscore[,1], col=3, lty=2)
lines(c(30,40,50,60,70), closedmeanscore[,3], col=3, lty=2)

lines(c(30,40,50,60,70), kohmeanscore[,2], lwd=2, col=7) # Yellow; KOH
lines(c(30,40,50,60,70), kohmeanscore[,1], col=7, lty=2)
lines(c(30,40,50,60,70), kohmeanscore[,3], col=7, lty=2)

### Median score ###
plot(c(30,40,50,60,70), closedmedscore[,2], type="l", lwd=2, col=3,  # Green; closed
     ylim=c(min(c(closedmedscore, kohmedscore)), max(c(closedmedscore, kohmedscore)))) 
lines(c(30,40,50,60,70), closedmedscore[,1], col=3, lty=2)
lines(c(30,40,50,60,70), closedmedscore[,3], col=3, lty=2)

lines(c(30,40,50,60,70), kohmedscore[,2], lwd=2, col=7) # Yellow; KOH
lines(c(30,40,50,60,70), kohmedscore[,1], col=7, lty=2)
lines(c(30,40,50,60,70), kohmedscore[,3], col=7, lty=2)

### Mean CRPS ###
plot(c(30,40,50,60,70), closedmeancrps[,2], type="l", lwd=2, col=3,  # Green; closed
     ylim=c(min(c(closedmeancrps, kohmeancrps)), max(c(closedmeancrps, kohmeancrps)))) 
lines(c(30,40,50,60,70), closedmeancrps[,1], col=3, lty=2)
lines(c(30,40,50,60,70), closedmeancrps[,3], col=3, lty=2)

lines(c(30,40,50,60,70), kohmeancrps[,2], lwd=2, col=7) # Yellow; KOH
lines(c(30,40,50,60,70), kohmeancrps[,1], col=7, lty=2)
lines(c(30,40,50,60,70), kohmeancrps[,3], col=7, lty=2)

### Median CRPS ###
plot(c(30,40,50,60,70), closedmedcrps[,2], type="l", lwd=2, col=3,  # Green; closed
     ylim=c(min(c(closedmedcrps, kohmedcrps)), max(c(closedmedcrps, kohmedcrps)))) 
lines(c(30,40,50,60,70), closedmedcrps[,1], col=3, lty=2)
lines(c(30,40,50,60,70), closedmedcrps[,3], col=3, lty=2)

lines(c(30,40,50,60,70), kohmedcrps[,2], lwd=2, col=7) # Yellow; KOH
lines(c(30,40,50,60,70), kohmedcrps[,1], col=7, lty=2)
lines(c(30,40,50,60,70), kohmedcrps[,3], col=7, lty=2)



