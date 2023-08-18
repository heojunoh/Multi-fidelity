library(lhs)
library(laGP)
library(plgp)
library(MuFiCokriging)

crps <- function(x, mu, sig2){ # The smaller, the better (0 to infinity)
  if(any(sig2==0)) sig2[sig2==0] <- eps
  -sqrt(sig2)*(1/sqrt(pi)-2*dnorm((x-mu)/sqrt(sig2))-(x-mu)/sqrt(sig2)*(2*pnorm((x-mu)/sqrt(sig2))-1))
}

### A2 function ###
f1 <- function(x)
{
  exp(-1.4*x) * cos(3.5*pi*x)
}

f2 <- function(x)
{ 
  y1 <- f1(x)
  
  log(y1^2+sqrt(x+3))
}

f3 <- function(x)
{ 
  y2 <- f2(x)
  
  -exp(y2)*sin(y2)
}

x <- seq(0,1,0.01)
f1(x)
f2(x)
f3(x)

plot(x, f1(x), type="l", lwd=2, col="red", ylim=c(min(c(f1(x), f2(x), f3(x))), max(c(f1(x), f2(x), f3(x)))))
lines(x, f2(x), col="orange",lwd=2,lty=1)
curve(f3(x),add=TRUE, col="green",lwd=2,lty=1) # high fidelity(TRUE); Black

### training data ###
n1 <- 15; n2 <- 10; n3 <- 5

rep <- 100
result.A2.rmse <- matrix(NA, rep, 2)
# result.A2.meanscore <- matrix(NA, rep, 5)
# result.A2.medscore <- matrix(NA, rep, 5)
result.A2.meancrps <- matrix(NA, rep, 2)
# result.A2.medcrps <- matrix(NA, rep, 5)
result.A2.comptime <- matrix(NA, rep, 2)
colnames(result.A2.rmse) <- c("closed", "Cokriging")
# colnames(result.A2.meanscore) <- c("single", "closed", "direct", "Cokriging", "KOH") # The larger, the better
# colnames(result.A2.medscore) <- c("single", "closed", "direct", "Cokriging", "KOH") # The larger, the better
colnames(result.A2.meancrps) <- c("closed", "Cokriging") # The smaller, the better
# colnames(result.A2.medcrps) <- c("single", "closed", "direct", "Cokriging", "KOH") # The smaller, the better
colnames(result.A2.comptime) <- c("closed", "Cokriging") # The smaller, the better


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
  
  ### test data ###
  x <- seq(0,1,0.01)
  
  
  ### closed ###
  tic.closed <- proc.time()[3]
  fit.closed <- RNAmf2(X1, y1, X2, y2, X3, y3, kernel="sqex", constant=TRUE)
  predy <- predRNAmf2(fit.closed, x)$mu
  predsig2 <- predRNAmf2(fit.closed, x)$sig2
  toc.closed <- proc.time()[3]
  
  
  ### Cokriging ###
  tic.cokm <- proc.time()[3]
  fit.muficokm <- MuFicokm(formula = list(~1,~1,~1), MuFidesign = NestDesign, covtype="gauss",
                           lower=eps, upper=0.1,
                           # coef.trend = list(0,c(0,0),c(0,0)),
                           response = list(y1,y2,y3), nlevel = 3)
  pred.muficokm <- predict(fit.muficokm, x, "SK")
  toc.cokm <- proc.time()[3]
  
  ### RMSE ###
  # result.A2.rmse[i,1] <- sqrt(sum((pred3$mu-f3(x))^2))/(sqrt(sum((f3(x))^2))) # single fidelity
  result.A2.rmse[i,1] <- sqrt(mean((predy-f3(x))^2)) # closed form
  # result.A2.rmse[i,3] <- sqrt(sum((pred3new$mu-f3(x))^2))/(sqrt(sum((f3(x))^2))) # not closed form
  result.A2.rmse[i,2] <- sqrt(mean((pred.muficokm$mean-f3(x))^2)) # Cokriging
  # result.A2.rmse[i,5] <- sqrt(sum((mx2-f3(x))^2))/(sqrt(sum((f3(x))^2))) # KOH

  # result.A2.meanscore[i,1] <- mean(score(f3(x), pred3$mu, pred3$sig2)) # single fidelity
  # result.A2.meanscore[i,2] <- mean(score(f3(x), predy, predsig2)) # closed form
  # result.A2.meanscore[i,3] <- mean(score(f3(x), pred3new$mu, pred3new$sig2)) # not closed form
  # result.A2.meanscore[i,4] <- mean(score(f3(x), pred.muficokm$mean, pred.muficokm$sig2)) # Cokriging
  # result.A2.meanscore[i,5] <- mean(score(f3(x), mx2, koh.var2)) # KOH
  # 
  # result.A2.medscore[i,1] <- median(score(f3(x), pred3$mu, pred3$sig2)) # single fidelity
  # result.A2.medscore[i,2] <- median(score(f3(x), predy, predsig2)) # closed form
  # result.A2.medscore[i,3] <- median(score(f3(x), pred3new$mu, pred3new$sig2)) # not closed form
  # result.A2.medscore[i,4] <- median(score(f3(x), pred.muficokm$mean, pred.muficokm$sig2)) # Cokriging
  # result.A2.medscore[i,5] <- median(score(f3(x), mx2, koh.var2)) # KOH
  
  # result.A2.meancrps[i,1] <- mean(crps(f3(x), pred3$mu, pred3$sig2)) # single fidelity
  result.A2.meancrps[i,1] <- mean(crps(f3(x), predy, predsig2)) # closed form
  # result.A2.meancrps[i,3] <- mean(crps(f3(x), pred3new$mu, pred3new$sig2)) # not closed form
  result.A2.meancrps[i,2] <- mean(crps(f3(x), pred.muficokm$mean, pred.muficokm$sig2)) # Cokriging
  # result.A2.meancrps[i,5] <- mean(crps(f3(x), mx2, koh.var2)) # KOH
  
  # result.A2.medcrps[i,1] <- median(crps(f3(x), pred3$mu, pred3$sig2)) # single fidelity
  # result.A2.medcrps[i,2] <- median(crps(f3(x), predy, predsig2)) # closed form
  # result.A2.medcrps[i,3] <- median(crps(f3(x), pred3new$mu, pred3new$sig2)) # not closed form
  # result.A2.medcrps[i,4] <- median(crps(f3(x), pred.muficokm$mean, pred.muficokm$sig2)) # Cokriging
  # result.A2.medcrps[i,5] <- median(crps(f3(x), mx2, koh.var2)) # KOH
  
  result.A2.comptime[i,1] <- toc.closed - tic.closed
  result.A2.comptime[i,2] <- toc.cokm - tic.cokm
}

par(mfrow=c(1,1))
#RMSE comparison#
apply(result.A2.rmse, 2, mean)
table(apply(result.A2.rmse, 1, which.min))
boxplot(result.A2.rmse)

# #score comparison, The larger, the better
# apply(result.A2.meanscore, 2, mean)
# table(apply(result.A2.meanscore, 1, which.max))
# boxplot(result.A2.meanscore)
# 
# #score comparison, The larger, the better
# apply(result.A2.medscore, 2, mean)
# table(apply(result.A2.medscore, 1, which.max))
# boxplot(result.A2.medscore)

#CRPS comparison, The smaller, the better
apply(result.A2.meancrps, 2, mean)
table(apply(result.A2.meancrps, 1, which.min))
boxplot(result.A2.meancrps)



