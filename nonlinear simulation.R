### nonlinear Example ###
library(lhs)
library(laGP)
source("GP.R")
source("KOH.R")
source("closed.R")
source("score.R")


fit.KOH <- function(X1, X2, Y1, Y2, g=eps){ # need to change function for another example
  
  ### KOH method ###
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
  
  return(list(b=c(b1, b2), rho=rho1, tau2hat=c(sig2_1, sig2_2), g=g, X1=X1, X2=X2, Y1=Y1, Y2=Y2))
}

pred.KOH <- function(fit, x){ # need to change function for another example
  
  X1 <- fit$X1
  X2 <- fit$X2
  Y1 <- fit$Y1
  Y2 <- fit$Y2
  
  b <- fit$b
  rho <- fit$rho
  tau2hat <- fit$tau2hat
  g <- fit$g
  
  
  ### prediction of 2nd order KOH ###
  tx1 <- cbind(rho*tau2hat[1]*covar.sep(x, X1, d=1/b[1], g=g), 
               rho^2*tau2hat[1]*covar.sep(x, X2, d=1/b[1], g=g) + tau2hat[2]*covar.sep(x, X2, d=1/b[2], g=g))
  
  V1 <- tau2hat[1]*covar.sep(X1, d=1/b[1], g=g)
  V12 <- rho*tau2hat[1]*covar.sep(X1, X2, d=1/b[1], g=0)
  V2 <- rho^2*tau2hat[1]*covar.sep(X2, d=1/b[1], g=g) + tau2hat[2]*covar.sep(X2, d=1/b[2], g=g)
  
  V_2 <- rbind(cbind(V1, V12), cbind(t(V12), V2))
  
  mx1 <- tx1 %*% solve(V_2) %*% c(Y1, Y2)
  
  ### posterior variance ###
  koh.var1 <- pmax(0, diag(tau2hat[2]*covar.sep(as.matrix(x), d=1/b[2], g=g) + tau2hat[1]*rho^2*covar.sep(as.matrix(x), d=1/b[1], g=g) - tx1 %*% solve(V_2)%*%t(tx1)))
  
  return(list(mu=mx1, sig2=koh.var1))
}

### synthetic function ###
f1 <- function(x)
{
  sin(8*pi*x)
}

f2 <- function(x)
{ 
  (x-sqrt(2))*(sin(8*pi*x))^2
}

x <- seq(0,1,length.out=400)
f1(x)
f2(x)

plot(x, f1(x), type="l", lwd=2, col="red", ylim=c(min(c(f1(x), f2(x))), max(c(f1(x), f2(x)))))
curve(f2(x),add=TRUE, col="green",lwd=2,lty=1) # high fidelity(TRUE); Black

### training data ###
n1 <- 35; n2 <- 30

rep <- 100
result.nonlinear.rmse <- matrix(NA, rep, 5)
result.nonlinear.meanscore <- matrix(NA, rep, 5)
result.nonlinear.medscore <- matrix(NA, rep, 5)
result.nonlinear.meancrps <- matrix(NA, rep, 5)
result.nonlinear.medcrps <- matrix(NA, rep, 5)
colnames(result.nonlinear.rmse) <- c("single", "closed", "direct", "Cokriging", "KOH")
colnames(result.nonlinear.meanscore) <- c("single", "closed", "direct", "Cokriging", "KOH") # The larger, the better
colnames(result.nonlinear.medscore) <- c("single", "closed", "direct", "Cokriging", "KOH") # The larger, the better
colnames(result.nonlinear.meancrps) <- c("single", "closed", "direct", "Cokriging", "KOH") # The smaller, the better
colnames(result.nonlinear.medcrps) <- c("single", "closed", "direct", "Cokriging", "KOH") # The smaller, the better

for(i in 1:rep) {
  set.seed(i)
  
  X1 <- maximinLHS(n1, 1)
  X2 <- maximinLHS(n2, 1)
  
  NestDesign <- NestedDesignBuild(design = list(X1,X2))
  
  X1 <- NestDesign$PX
  X2 <- ExtractNestDesign(NestDesign,2)
  
  y1 <- f1(X1)
  y2 <- f2(X2)
  
  
  ### model fitting for f1 ###
  eps <- sqrt(.Machine$double.eps)
  fit.GP1 <- GP(X1, y1, constant=TRUE)
  
  ### model fitting using (x2, f1(x2)) ###
  w1.x2 <- pred.GP(fit.GP1, X2)$mu # can interpolate; nested
  X2new <- cbind(X2, w1.x2) # combine (X2, f1(x2))
  fit.GP2new <- GP(X2new, y2, constant=TRUE) # model fitting for f_M(X2, f1(x2))
  
  ### test data ###
  x <- seq(0,1,0.01)
  
  
  ### closed ###
  predy <- closed(x, fit.GP1, fit.GP2new, constant=TRUE)$mu
  predsig2 <- closed(x, fit.GP1, fit.GP2new, constant=TRUE)$sig2
  
  
  # ### KOH method ###
  # fit.KOH2 <- fit.KOH(X1, X2, y1, y2)
  # pred.KOH2 <- pred.KOH(fit.KOH2, x)
  # mx1 <- pred.KOH2$mu
  # koh.var1 <- pred.KOH2$sig2
  # 
  # 
  # ### direct fitting; not using closed form. f1(u) and f_M(u) from (u, f_M(u, f1(u))) are random variables.
  # # w1.x <- rnorm(length(x), mean=pred.GP(fit.GP1, x)$mu, sd=sqrt(pred.GP(fit.GP1, x)$sig2)) # sample f1(x)
  # # w2.x <- rnorm(length(x), mean=pred.GP(fit.GP2new, cbind(x, w1.x))$mu, sd=sqrt(pred.GP(fit.GP2new, cbind(x, w1.x))$sig2))
  # # xxnew <- cbind(x, w2.x)
  # # pred3new <- pred.GP(fit.GP3new, xxnew) # not closed form
  # w1.x <- c(rep(NA, length(x)))
  # w2.x <- c(rep(NA, length(x)))
  # for(j in 1:length(x)){
  #   w1.x[j] <- mean(rnorm(10000, mean=pred.GP(fit.GP1, x[j])$mu, sd=sqrt(pred.GP(fit.GP1, x[j])$sig2)))
  # }
  # 
  # xxnew <- cbind(x, w1.x)
  # pred2new <- pred.GP(fit.GP2new, xxnew) # not closed form
  # 
  # ### prediction of original GP with single fidelity ###
  # fit.GP2 <- GP(X2, y2, constant=TRUE)
  # pred2 <- pred.GP(fit.GP2, x)
  # 
  # ### Cokriging ###
  # fit.muficokm <- MuFicokm(formula = list(~1,~1), MuFidesign = NestDesign, #covtype="gauss",
  #                          coef.trend = list(0,c(0,0)), response = list(y1,y2), nlevel = 2)
  # pred.muficokm <- predict(fit.muficokm, x, "SK")
  
  ### RMSE ###
  result.nonlinear.rmse[i,1] <- sqrt(sum((pred2$mu-f2(x))^2))/(sqrt(sum((f2(x))^2))) # single fidelity
  result.nonlinear.rmse[i,2] <- sqrt(sum((predy-f2(x))^2))/(sqrt(sum((f2(x))^2))) # closed form
  # result.nonlinear.rmse[i,3] <- sqrt(sum((pred2new$mu-f2(x))^2))/(sqrt(sum((f2(x))^2))) # not closed form
  # result.nonlinear.rmse[i,4] <- sqrt(sum((pred.muficokm$mean-f2(x))^2))/(sqrt(sum((f2(x))^2))) # Cokriging
  # result.nonlinear.rmse[i,5] <- sqrt(sum((mx1-f2(x))^2))/(sqrt(sum((f2(x))^2))) # KOH
  # 
  result.nonlinear.meanscore[i,1] <- mean(score(f2(x), pred2$mu, pred2$sig2)) # single fidelity
  result.nonlinear.meanscore[i,2] <- mean(score(f2(x), predy, predsig2)) # closed form
  # result.nonlinear.meanscore[i,3] <- mean(score(f2(x), pred2new$mu, pred2new$sig2)) # not closed form
  # result.nonlinear.meanscore[i,4] <- mean(score(f2(x), pred.muficokm$mean, pred.muficokm$sig2)) # Cokriging
  # result.nonlinear.meanscore[i,5] <- mean(score(f2(x), mx1, koh.var1)) # KOH
  # 
  result.nonlinear.medscore[i,1] <- median(score(f2(x), pred2$mu, pred2$sig2)) # single fidelity
  result.nonlinear.medscore[i,2] <- median(score(f2(x), predy, predsig2)) # closed form
  # result.nonlinear.medscore[i,3] <- median(score(f2(x), pred2new$mu, pred2new$sig2)) # not closed form
  # result.nonlinear.medscore[i,4] <- median(score(f2(x), pred.muficokm$mean, pred.muficokm$sig2)) # Cokriging
  # result.nonlinear.medscore[i,5] <- median(score(f2(x), mx1, koh.var1)) # KOH
  
  result.nonlinear.meancrps[i,1] <- mean(crps(f2(x), pred2$mu, pred2$sig2)) # single fidelity
  result.nonlinear.meancrps[i,2] <- mean(crps(f2(x), predy, predsig2)) # closed form
  # result.nonlinear.meancrps[i,3] <- mean(crps(f2(x), pred2new$mu, pred2new$sig2)) # not closed form
  # result.nonlinear.meancrps[i,4] <- mean(crps(f2(x), pred.muficokm$mean, pred.muficokm$sig2)) # Cokriging
  # result.nonlinear.meancrps[i,5] <- mean(crps(f2(x), mx1, koh.var1)) # KOH
  
  result.nonlinear.medcrps[i,1] <- median(crps(f2(x), pred2$mu, pred2$sig2)) # single fidelity
  result.nonlinear.medcrps[i,2] <- median(crps(f2(x), predy, predsig2)) # closed form
  # result.nonlinear.medcrps[i,3] <- median(crps(f2(x), pred2new$mu, pred2new$sig2)) # not closed form
  # result.nonlinear.medcrps[i,4] <- median(crps(f2(x), pred.muficokm$mean, pred.muficokm$sig2)) # Cokriging
  # result.nonlinear.medcrps[i,5] <- median(crps(f2(x), mx1, koh.var1)) # KOH
  
}


par(mfrow=c(1,1))
#RMSE comparison#
apply(result.nonlinear.rmse, 2, mean)
table(apply(result.nonlinear.rmse, 1, which.min))
boxplot(result.nonlinear.rmse)
log(apply(result.nonlinear.rmse, 2, mean))
#score comparison, The larger, the better
apply(result.nonlinear.meanscore, 2, mean)
table(apply(result.nonlinear.meanscore, 1, which.max))
boxplot(result.nonlinear.meanscore)

#score comparison, The larger, the better
apply(result.nonlinear.medscore, 2, mean)
table(apply(result.nonlinear.medscore, 1, which.max))
boxplot(result.nonlinear.medscore)

#CRPS comparison, The smaller, the better
apply(result.nonlinear.meancrps, 2, mean)
table(apply(result.nonlinear.meancrps, 1, which.min))
boxplot(result.nonlinear.meancrps)

#CRPS comparison, The smaller, the better
apply(result.nonlinear.medcrps, 2, mean)
table(apply(result.nonlinear.medcrps, 1, which.min))
boxplot(result.nonlinear.medcrps)


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


closedrmse <- rbind(closedrmse, c(sort(result.nonlinear.rmse[,2])[5], apply(result.nonlinear.rmse, 2, mean)[2], sort(result.nonlinear.rmse[,2])[96]))
closedmeanscore <- rbind(closedmeanscore, c(sort(result.nonlinear.meanscore[,2])[5], apply(result.nonlinear.meanscore, 2, mean)[2], sort(result.nonlinear.meanscore[,2])[96]))
closedmedscore <- rbind(closedmedscore, c(sort(result.nonlinear.medscore[,2])[5], apply(result.nonlinear.medscore, 2, mean)[2], sort(result.nonlinear.medscore[,2])[96]))
closedmeancrps <- rbind(closedmeancrps, c(sort(result.nonlinear.meancrps[,2])[5], apply(result.nonlinear.meancrps, 2, mean)[2], sort(result.nonlinear.meancrps[,2])[96]))
closedmedcrps <- rbind(closedmedcrps, c(sort(result.nonlinear.medcrps[,2])[5], apply(result.nonlinear.medcrps, 2, mean)[2], sort(result.nonlinear.medcrps[,2])[96]))

kohrmse <- rbind(kohrmse, c(sort(result.nonlinear.rmse[,4])[5], apply(result.nonlinear.rmse, 2, mean)[4], sort(result.nonlinear.rmse[,4])[96]))
kohmeanscore <- rbind(kohmeanscore, c(sort(result.nonlinear.meanscore[,4])[5], apply(result.nonlinear.meanscore, 2, mean)[4], sort(result.nonlinear.meanscore[,4])[96]))
kohmedscore <- rbind(kohmedscore, c(sort(result.nonlinear.medscore[,4])[5], apply(result.nonlinear.medscore, 2, mean)[4], sort(result.nonlinear.medscore[,4])[96]))
kohmeancrps <- rbind(kohmeancrps, c(sort(result.nonlinear.meancrps[,4])[5], apply(result.nonlinear.meancrps, 2, mean)[4], sort(result.nonlinear.meancrps[,4])[96]))
kohmedcrps <- rbind(kohmedcrps, c(sort(result.nonlinear.medcrps[,4])[5], apply(result.nonlinear.medcrps, 2, mean)[4], sort(result.nonlinear.medcrps[,4])[96]))


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



