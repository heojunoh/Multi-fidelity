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
  Y2d3 <- apply(X3,1,output.g1, l=3)
  Y1d2 <- apply(X2,1,output.g1, l=1)
  
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
Gramacy1 <- function(xx, l){
  x1 <- xx[1]
  x2 <- xx[2]
  
  # 10*x1*exp(-x1^2-x2^2)
  10*x1*exp(-x1^2-x2^2)* (1+0.8^l) 
  }

output.g1 <- function(x, l){
  factor_range <- list("x1" = c(-2, 4), "x2" = c(-2, 4))
  
  for(i in 1:length(factor_range)) x[i] <- factor_range[[i]][1] + x[i] * diff(factor_range[[i]])
  Gramacy1(x[1:2], l)
} 

### training data ###
n1 <- 15; n2 <- 11; n3 <- 7

rep <- 100
result.gramacy1.rmse <- matrix(NA, rep, 4)
result.gramacy1.meanscore <- matrix(NA, rep, 4)
result.gramacy1.medscore <- matrix(NA, rep, 4)
result.gramacy1.meancrps <- matrix(NA, rep, 4)
result.gramacy1.medcrps <- matrix(NA, rep, 4)
colnames(result.gramacy1.rmse) <- c("single", "closed", "direct", "KOH")
colnames(result.gramacy1.meanscore) <- c("single", "closed", "direct", "KOH") # The larger, the better
colnames(result.gramacy1.medscore) <- c("single", "closed", "direct", "KOH") # The larger, the better
colnames(result.gramacy1.meancrps) <- c("single", "closed", "direct", "KOH") # The smaller, the better
colnames(result.gramacy1.medcrps) <- c("single", "closed", "direct", "KOH") # The smaller, the better


for(i in 1:rep) {
  set.seed(i)
  
  X1 <- maximinLHS(n1, 2)
  X2 <- maximinLHS(n2, 2)
  X3 <- maximinLHS(n3, 2)
  
  NestDesign <- NestedDesignBuild(design = list(X1,X2,X3))
  
  X1 <- NestDesign$PX
  X2 <- ExtractNestDesign(NestDesign,2)
  X3 <- ExtractNestDesign(NestDesign,3)
  
  y3 <- apply(X3,1,output.g1, l=5)
  y2 <- apply(X2,1,output.g1, l=3)
  y1 <- apply(X1,1,output.g1, l=1)
  
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
  
  # ### Cokriging ###
  # fit.muficokm <- MuFicokm(formula = list(~1,~1,~1), MuFidesign = NestDesign, #covtype="gauss",
  #                          coef.trend = list(0,c(0,0)), response = list(y1,y2,y3), nlevel = 3)
  # pred.muficokm <- predict(fit.muficokm, x, "UK")
  
  ### RMSE ###
  result.gramacy1.rmse[i,1] <- sqrt(mean((pred3$mu-apply(x,1,output.g1, l=Inf))^2)) # single fidelity
  result.gramacy1.rmse[i,2] <- sqrt(mean((predy-apply(x,1,output.g1, l=Inf))^2)) # closed form
  result.gramacy1.rmse[i,3] <- sqrt(mean((pred3new$mu-apply(x,1,output.g1, l=Inf))^2)) # not closed form
  # result.gramacy1.rmse[i,4] <- sqrt(mean((pred.muficokm$mean-apply(x,1,output.g1, l=Inf))^2)) # Cokriging
  result.gramacy1.rmse[i,4] <- sqrt(mean((mx2-apply(x,1,output.g1, l=Inf))^2)) # KOH
  
  result.gramacy1.meanscore[i,1] <- mean(score(apply(x,1,output.g1, l=Inf), pred3$mu, pred3$sig2)) # single fidelity
  result.gramacy1.meanscore[i,2] <- mean(score(apply(x,1,output.g1, l=Inf), predy, predsig2)) # closed form
  result.gramacy1.meanscore[i,3] <- mean(score(apply(x,1,output.g1, l=Inf), pred3new$mu, pred3new$sig2)) # not closed form
  # result.gramacy1.meanscore[i,4] <- mean(score(apply(x,1,output.g1, l=Inf), pred.muficokm$mean, pred.muficokm$sig2)) # Cokriging
  result.gramacy1.meanscore[i,4] <- mean(score(apply(x,1,output.g1, l=Inf), mx2, koh.var2)) # KOH
  
  result.gramacy1.medscore[i,1] <- median(score(apply(x,1,output.g1, l=Inf), pred3$mu, pred3$sig2)) # single fidelity
  result.gramacy1.medscore[i,2] <- median(score(apply(x,1,output.g1, l=Inf), predy, predsig2)) # closed form
  result.gramacy1.medscore[i,3] <- median(score(apply(x,1,output.g1, l=Inf), pred3new$mu, pred3new$sig2)) # not closed form
  # result.gramacy1.medscore[i,4] <- median(score(apply(x,1,output.g1, l=Inf), pred.muficokm$mean, pred.muficokm$sig2)) # Cokriging
  result.gramacy1.medscore[i,4] <- median(score(apply(x,1,output.g1, l=Inf), mx2, koh.var2)) # KOH
  
  result.gramacy1.meancrps[i,1] <- mean(crps(apply(x,1,output.g1, l=Inf), pred3$mu, pred3$sig2)) # single fidelity
  result.gramacy1.meancrps[i,2] <- mean(crps(apply(x,1,output.g1, l=Inf), predy, predsig2)) # closed form
  result.gramacy1.meancrps[i,3] <- mean(crps(apply(x,1,output.g1, l=Inf), pred3new$mu, pred3new$sig2)) # not closed form
  # result.gramacy1.meancrps[i,4] <- mean(crps(apply(x,1,output.g1, l=Inf), pred.muficokm$mean, pred.muficokm$sig2)) # Cokriging
  result.gramacy1.meancrps[i,4] <- mean(crps(apply(x,1,output.g1, l=Inf), mx2, koh.var2)) # KOH
  
  result.gramacy1.medcrps[i,1] <- median(crps(apply(x,1,output.g1, l=Inf), pred3$mu, pred3$sig2)) # single fidelity
  result.gramacy1.medcrps[i,2] <- median(crps(apply(x,1,output.g1, l=Inf), predy, predsig2)) # closed form
  result.gramacy1.medcrps[i,3] <- median(crps(apply(x,1,output.g1, l=Inf), pred3new$mu, pred3new$sig2)) # not closed form
  # result.gramacy1.medcrps[i,4] <- median(crps(apply(x,1,output.g1, l=Inf), pred.muficokm$mean, pred.muficokm$sig2)) # Cokriging
  result.gramacy1.medcrps[i,4] <- median(crps(apply(x,1,output.g1, l=Inf), mx2, koh.var2)) # KOH
  
}

par(mfrow=c(1,1))
#RMSE comparison#
apply(result.gramacy1.rmse, 2, mean)
table(apply(result.gramacy1.rmse, 1, which.min))
boxplot(result.gramacy1.rmse)

#score comparison, The larger, the better
apply(result.gramacy1.meanscore, 2, mean)
table(apply(result.gramacy1.meanscore, 1, which.max))
boxplot(result.gramacy1.meanscore)

#score comparison, The larger, the better
apply(result.gramacy1.medscore, 2, mean)
table(apply(result.gramacy1.medscore, 1, which.max))
boxplot(result.gramacy1.medscore)

#CRPS comparison, The smaller, the better
apply(result.gramacy1.meancrps, 2, mean)
table(apply(result.gramacy1.meancrps, 1, which.min))
boxplot(result.gramacy1.meancrps)

#CRPS comparison, The smaller, the better
apply(result.gramacy1.medcrps, 2, mean)
table(apply(result.gramacy1.medcrps, 1, which.min))
boxplot(result.gramacy1.medcrps)


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


closedrmse <- rbind(closedrmse, c(sort(result.gramacy1.rmse[,2])[5], apply(result.gramacy1.rmse, 2, mean)[2], sort(result.gramacy1.rmse[,2])[96]))
closedmeanscore <- rbind(closedmeanscore, c(sort(result.gramacy1.meanscore[,2])[5], apply(result.gramacy1.meanscore, 2, mean)[2], sort(result.gramacy1.meanscore[,2])[96]))
closedmedscore <- rbind(closedmedscore, c(sort(result.gramacy1.medscore[,2])[5], apply(result.gramacy1.medscore, 2, mean)[2], sort(result.gramacy1.medscore[,2])[96]))
closedmeancrps <- rbind(closedmeancrps, c(sort(result.gramacy1.meancrps[,2])[5], apply(result.gramacy1.meancrps, 2, mean)[2], sort(result.gramacy1.meancrps[,2])[96]))
closedmedcrps <- rbind(closedmedcrps, c(sort(result.gramacy1.medcrps[,2])[5], apply(result.gramacy1.medcrps, 2, mean)[2], sort(result.gramacy1.medcrps[,2])[96]))

kohrmse <- rbind(kohrmse, c(sort(result.gramacy1.rmse[,4])[5], apply(result.gramacy1.rmse, 2, mean)[4], sort(result.gramacy1.rmse[,4])[96]))
kohmeanscore <- rbind(kohmeanscore, c(sort(result.gramacy1.meanscore[,4])[5], apply(result.gramacy1.meanscore, 2, mean)[4], sort(result.gramacy1.meanscore[,4])[96]))
kohmedscore <- rbind(kohmedscore, c(sort(result.gramacy1.medscore[,4])[5], apply(result.gramacy1.medscore, 2, mean)[4], sort(result.gramacy1.medscore[,4])[96]))
kohmeancrps <- rbind(kohmeancrps, c(sort(result.gramacy1.meancrps[,4])[5], apply(result.gramacy1.meancrps, 2, mean)[4], sort(result.gramacy1.meancrps[,4])[96]))
kohmedcrps <- rbind(kohmedcrps, c(sort(result.gramacy1.medcrps[,4])[5], apply(result.gramacy1.medcrps, 2, mean)[4], sort(result.gramacy1.medcrps[,4])[96]))


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
plot(c(7:20), closedrmse[,2], type="l", lwd=2, col=3,  # Green; closed
     ylim=c(min(c(closedrmse, kohrmse)), max(c(closedrmse, kohrmse)))) 
lines(c(7:20), closedrmse[,1], col=3, lty=2)
lines(c(7:20), closedrmse[,3], col=3, lty=2)

lines(c(7:20), kohrmse[,2], lwd=2, col=7) # Yellow; KOH
lines(c(7:20), kohrmse[,1], col=7, lty=2)
lines(c(7:20), kohrmse[,3], col=7, lty=2)

### Mean score ###
plot(c(7:20), closedmeanscore[,2], type="l", lwd=2, col=3,  # Green; closed
     ylim=c(min(c(closedmeanscore, kohmeanscore)), max(c(closedmeanscore, kohmeanscore)))) 
lines(c(7:20), closedmeanscore[,1], col=3, lty=2)
lines(c(7:20), closedmeanscore[,3], col=3, lty=2)

lines(c(7:20), kohmeanscore[,2], lwd=2, col=7) # Yellow; KOH
lines(c(7:20), kohmeanscore[,1], col=7, lty=2)
lines(c(7:20), kohmeanscore[,3], col=7, lty=2)

### Median score ###
plot(c(7:20), closedmedscore[,2], type="l", lwd=2, col=3,  # Green; closed
     ylim=c(min(c(closedmedscore, kohmedscore)), max(c(closedmedscore, kohmedscore)))) 
lines(c(7:20), closedmedscore[,1], col=3, lty=2)
lines(c(7:20), closedmedscore[,3], col=3, lty=2)

lines(c(7:20), kohmedscore[,2], lwd=2, col=7) # Yellow; KOH
lines(c(7:20), kohmedscore[,1], col=7, lty=2)
lines(c(7:20), kohmedscore[,3], col=7, lty=2)

### Mean CRPS ###
plot(c(7:20), closedmeancrps[,2], type="l", lwd=2, col=3,  # Green; closed
     ylim=c(min(c(closedmeancrps, kohmeancrps)), max(c(closedmeancrps, kohmeancrps)))) 
lines(c(7:20), closedmeancrps[,1], col=3, lty=2)
lines(c(7:20), closedmeancrps[,3], col=3, lty=2)

lines(c(7:20), kohmeancrps[,2], lwd=2, col=7) # Yellow; KOH
lines(c(7:20), kohmeancrps[,1], col=7, lty=2)
lines(c(7:20), kohmeancrps[,3], col=7, lty=2)

### Median CRPS ###
plot(c(7:20), closedmedcrps[,2], type="l", lwd=2, col=3,  # Green; closed
     ylim=c(min(c(closedmedcrps, kohmedcrps)), max(c(closedmedcrps, kohmedcrps)))) 
lines(c(7:20), closedmedcrps[,1], col=3, lty=2)
lines(c(7:20), closedmedcrps[,3], col=3, lty=2)

lines(c(7:20), kohmedcrps[,2], lwd=2, col=7) # Yellow; KOH
lines(c(7:20), kohmedcrps[,1], col=7, lty=2)
lines(c(7:20), kohmedcrps[,3], col=7, lty=2)


