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
  
  term1 + term2*5*0.8^l + (term1+term2)^3 + exp(-2*term1*term2)
}

### training data ###
n1 <- 11; n2 <- 9; n3 <- 7

rep <- 100
result.3layers.rmse <- matrix(NA, rep, 5)
result.3layers.meanscore <- matrix(NA, rep, 5)
result.3layers.medscore <- matrix(NA, rep, 5)
result.3layers.meancrps <- matrix(NA, rep, 5)
result.3layers.medcrps <- matrix(NA, rep, 5)
colnames(result.3layers.rmse) <- c("single", "closed", "direct", "Cokriging", "KOH")
colnames(result.3layers.meanscore) <- c("single", "closed", "direct", "Cokriging", "KOH") # The larger, the better
colnames(result.3layers.medscore) <- c("single", "closed", "direct", "Cokriging", "KOH") # The larger, the better
colnames(result.3layers.meancrps) <- c("single", "closed", "direct", "Cokriging", "KOH") # The smaller, the better
colnames(result.3layers.medcrps) <- c("single", "closed", "direct", "Cokriging", "KOH") # The smaller, the better


for(i in 1:rep) {
  set.seed(i)
  
  # X3 <- maximinLHS(n3, 1) # x^H
  # y3 <- fl(X3, l=5)
  # X2 <- matrix(c(X3, maximinLHS(n2-n3, 1))) # x^M
  # y2 <- fl(X2, l=3)
  # X1 <- matrix(c(X2, maximinLHS(n1-n2, 1))) # x^L
  # y1 <- fl(X1, l=1)
  
  X1 <- maximinLHS(n1, 1)
  X2 <- maximinLHS(n2, 1)
  X3 <- maximinLHS(n3, 1)
  
  NestDesign <- NestedDesignBuild(design = list(X1,X2,X3))
  
  X1 <- NestDesign$PX
  X2 <- ExtractNestDesign(NestDesign,2)
  X3 <- ExtractNestDesign(NestDesign,3)
  
  y1 <- fl(X1, l=1)
  y2 <- fl(X2, l=3)
  y3 <- fl(X3, l=5)
  
  ### model fitting for f1 ###
  eps <- sqrt(.Machine$double.eps)
  fit.GP1 <- GP(X1, y1, constant=TRUE)
  
  ### model fitting using (x2, f1(x2)) ###
  w1.x2 <- pred.GP(fit.GP1, X2)$mu # can interpolate; nested
  X2new <- cbind(X2, w1.x2) # combine (X2, f1(x2))
  fit.GP2new <- GP(X2new, y2, constant=TRUE) # model fitting for f_M(X2, f1(x2))
  
  ### model fitting using (x3, f2(x3, f1(x3))) ###
  w1.x3 <- pred.GP(fit.GP1, X3)$mu # can interpolate; nested
  w2.x3 <- pred.GP(fit.GP2new, cbind(X3, w1.x3))$mu # can interpolate; nested
  X3new <- cbind(X3, w2.x3) # combine (X3, f2(x3, f1(x3)))
  fit.GP3new <- GP(X3new, y3, constant=TRUE) # model fitting for f_H(X3, f2(x3, f1(x3)))
  
  
  ### test data ###
  x <- seq(0,1,0.01)
  
  
  ### closed ###
  predy <- closed2(x, fit.GP1, fit.GP2new, fit.GP3new, constant=TRUE)$mu
  predsig2 <- closed2(x, fit.GP1, fit.GP2new, fit.GP3new, constant=TRUE)$sig2
  
  
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
  fit.GP3 <- GP(X3, y3, constant=TRUE)
  pred3 <- pred.GP(fit.GP3, x)
  
  ### Cokriging ###
  fit.muficokm <- MuFicokm(formula = list(~1,~1,~1), MuFidesign = NestDesign, #covtype="gauss",
                           coef.trend = list(0,c(0,0),c(0,0)), response = list(y1,y2,y3), nlevel = 3)
  pred.muficokm <- predict(fit.muficokm, x, "SK")
  
  ### RMSE ###
  result.3layers.rmse[i,1] <- sqrt(mean((pred3$mu-fl(x, l=Inf))^2)) # single fidelity
  result.3layers.rmse[i,2] <- sqrt(mean((predy-fl(x, l=Inf))^2)) # closed form
  result.3layers.rmse[i,3] <- sqrt(mean((pred3new$mu-fl(x, l=Inf))^2)) # not closed form
  result.3layers.rmse[i,4] <- sqrt(mean((pred.muficokm$mean-fl(x, l=Inf))^2)) # Cokriging
  result.3layers.rmse[i,5] <- sqrt(mean((mx2-fl(x, l=Inf))^2)) # KOH
  
  result.3layers.meanscore[i,1] <- mean(score(fl(x, l=Inf), pred3$mu, pred3$sig2)) # single fidelity
  result.3layers.meanscore[i,2] <- mean(score(fl(x, l=Inf), predy, predsig2)) # closed form
  result.3layers.meanscore[i,3] <- mean(score(fl(x, l=Inf), pred3new$mu, pred3new$sig2)) # not closed form
  result.3layers.meanscore[i,4] <- mean(score(fl(x, l=Inf), pred.muficokm$mean, pred.muficokm$sig2)) # Cokriging
  result.3layers.meanscore[i,5] <- mean(score(fl(x, l=Inf), mx2, koh.var2)) # KOH
  
  result.3layers.medscore[i,1] <- median(score(fl(x, l=Inf), pred3$mu, pred3$sig2)) # single fidelity
  result.3layers.medscore[i,2] <- median(score(fl(x, l=Inf), predy, predsig2)) # closed form
  result.3layers.medscore[i,3] <- median(score(fl(x, l=Inf), pred3new$mu, pred3new$sig2)) # not closed form
  result.3layers.medscore[i,4] <- median(score(fl(x, l=Inf), pred.muficokm$mean, pred.muficokm$sig2)) # Cokriging
  result.3layers.medscore[i,5] <- median(score(fl(x, l=Inf), mx2, koh.var2)) # KOH
  
  result.3layers.meancrps[i,1] <- mean(crps(fl(x, l=Inf), pred3$mu, pred3$sig2)) # single fidelity
  result.3layers.meancrps[i,2] <- mean(crps(fl(x, l=Inf), predy, predsig2)) # closed form
  result.3layers.meancrps[i,3] <- mean(crps(fl(x, l=Inf), pred3new$mu, pred3new$sig2)) # not closed form
  result.3layers.meancrps[i,4] <- mean(crps(fl(x, l=Inf), pred.muficokm$mean, pred.muficokm$sig2)) # Cokriging
  result.3layers.meancrps[i,5] <- mean(crps(fl(x, l=Inf), mx2, koh.var2)) # KOH
  
  result.3layers.medcrps[i,1] <- median(crps(fl(x, l=Inf), pred3$mu, pred3$sig2)) # single fidelity
  result.3layers.medcrps[i,2] <- median(crps(fl(x, l=Inf), predy, predsig2)) # closed form
  result.3layers.medcrps[i,3] <- median(crps(fl(x, l=Inf), pred3new$mu, pred3new$sig2)) # not closed form
  result.3layers.medcrps[i,4] <- median(crps(fl(x, l=Inf), pred.muficokm$mean, pred.muficokm$sig2)) # Cokriging
  result.3layers.medcrps[i,5] <- median(crps(fl(x, l=Inf), mx2, koh.var2)) # KOH
  
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


closedrmse <- rbind(closedrmse, c(sort(result.3layers.rmse[,2])[5], apply(result.3layers.rmse, 2, mean)[2], sort(result.3layers.rmse[,2])[96]))
closedmeanscore <- rbind(closedmeanscore, c(sort(result.3layers.meanscore[,2])[5], apply(result.3layers.meanscore, 2, mean)[2], sort(result.3layers.meanscore[,2])[96]))
closedmedscore <- rbind(closedmedscore, c(sort(result.3layers.medscore[,2])[5], apply(result.3layers.medscore, 2, mean)[2], sort(result.3layers.medscore[,2])[96]))
closedmeancrps <- rbind(closedmeancrps, c(sort(result.3layers.meancrps[,2])[5], apply(result.3layers.meancrps, 2, mean)[2], sort(result.3layers.meancrps[,2])[96]))
closedmedcrps <- rbind(closedmedcrps, c(sort(result.3layers.medcrps[,2])[5], apply(result.3layers.medcrps, 2, mean)[2], sort(result.3layers.medcrps[,2])[96]))

kohrmse <- rbind(kohrmse, c(sort(result.3layers.rmse[,5])[5], apply(result.3layers.rmse, 2, mean)[5], sort(result.3layers.rmse[,5])[96]))
kohmeanscore <- rbind(kohmeanscore, c(sort(result.3layers.meanscore[,5])[5], apply(result.3layers.meanscore, 2, mean)[5], sort(result.3layers.meanscore[,5])[96]))
kohmedscore <- rbind(kohmedscore, c(sort(result.3layers.medscore[,5])[5], apply(result.3layers.medscore, 2, mean)[5], sort(result.3layers.medscore[,5])[96]))
kohmeancrps <- rbind(kohmeancrps, c(sort(result.3layers.meancrps[,5])[5], apply(result.3layers.meancrps, 2, mean)[5], sort(result.3layers.meancrps[,5])[96]))
kohmedcrps <- rbind(kohmedcrps, c(sort(result.3layers.medcrps[,5])[5], apply(result.3layers.medcrps, 2, mean)[5], sort(result.3layers.medcrps[,5])[96]))


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


