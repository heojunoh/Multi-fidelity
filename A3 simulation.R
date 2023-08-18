### A3 Example ###
library(lhs)
library(laGP)
library(MuFiCokriging)

crps <- function(x, mu, sig2){ # The smaller, the better (0 to infinity)
  if(any(sig2==0)) sig2[sig2==0] <- eps
  -sqrt(sig2)*(1/sqrt(pi)-2*dnorm((x-mu)/sqrt(sig2))-(x-mu)/sqrt(sig2)*(2*pnorm((x-mu)/sqrt(sig2))-1))
}


### synthetic function ###
A3low <- function(xx){
  x1 <- xx[1]
  x2 <- xx[2]
  
  exp(-1.4*(x1*exp(-x1^2-x2^2)+0.5))*cos(3.5*pi*(x1*exp(-x1^2-x2^2)+0.5))
}

A3high <- function(xx)
{ 
  x1 <- xx[1]
  x2 <- xx[2]
  
  x1*exp(-x1^2-x2^2)
}


output.A3low <- function(x){
  factor_range <- list("x1" = c(-2, 6), "x2" = c(-2, 6))
  
  for(i in 1:length(factor_range)) x[i] <- factor_range[[i]][1] + x[i] * diff(factor_range[[i]])
  A3low(x[1:2])
} 

output.A3high <- function(x){
  factor_range <- list("x1" = c(-2, 6), "x2" = c(-2, 6))
  
  for(i in 1:length(factor_range)) x[i] <- factor_range[[i]][1] + x[i] * diff(factor_range[[i]])
  A3high(x[1:2])
} 



### training data ###
n1 <- 20; n2 <- 10

rep <- 100
result.A3.rmse <- matrix(NA, rep, 2)
result.A3.meancrps <- matrix(NA, rep, 2)
result.A3.comptime <- matrix(NA, rep, 2)
colnames(result.A3.rmse) <- c("closed", "Cokriging")
colnames(result.A3.meancrps) <- c("closed", "Cokriging") # The smaller, the better
colnames(result.A3.comptime) <- c("closed", "Cokriging") # The smaller, the better
meanf <- c(rep(0,100))

for(i in 1:rep) {
  set.seed(i)
  
  X1 <- maximinLHS(n1, 2)
  X2 <- maximinLHS(n2, 2)
  
  NestDesign <- NestedDesignBuild(design = list(X1,X2))
  
  X1 <- NestDesign$PX
  X2 <- ExtractNestDesign(NestDesign,2)
  
  y1 <- apply(X1,1,output.A3low)
  y2 <- apply(X2,1,output.A3high)
  
  
  # ### model fitting for f1 ###
  eps <- sqrt(.Machine$double.eps)
  # fit.GP1 <- GP(X1, y1, constant=TRUE)
  # 
  # ### model fitting using (x2, f1(x2)) ###
  # w1.x2 <- pred.GP(fit.GP1, X2)$mu # can interpolate; nested
  # X2new <- cbind(X2, w1.x2) # combine (X2, f1(x2))
  # fit.GP2new <- GP(X2new, y2, constant=TRUE) # model fitting for f_M(X2, f1(x2))
  # 
  # ### model fitting using (x3, f2(x3, f1(x3))) ###
  # w1.x3 <- pred.GP(fit.GP1, X3)$mu # can interpolate; nested
  # w2.x3 <- pred.GP(fit.GP2new, cbind(X3, w1.x3))$mu # can interpolate; nested
  # X3new <- cbind(X3, w2.x3) # combine (X3, f2(x3, f1(x3)))
  # fit.GP3new <- GP(X3new, y3, constant=TRUE) # model fitting for f_H(X3, f2(x3, f1(x3)))
  
  
  ### test data ###
  x <- maximinLHS(100, 2)
  
  ### closed ###
  tic.closed <- proc.time()[3]
  fit.closed <- RNAmf(X1, y1, X2, y2, kernel="sqex", constant=TRUE)
  predy <- predRNAmf(fit.closed, x)$mu
  predsig2 <- predRNAmf(fit.closed, x)$sig2
  toc.closed <- proc.time()[3]

  
  ### Cokriging ###
  tic.cokm <- proc.time()[3]
  fit.muficokm <- MuFicokm(formula = list(~1,~1), MuFidesign = NestDesign, covtype="gauss",
                           lower=eps, upper=0.5,
                           # coef.trend = list(0,c(0,0),c(0,0)),
                           response = list(y1,y2), nlevel = 2)
  pred.muficokm <- predict(fit.muficokm, x, "SK")
  toc.cokm <- proc.time()[3]
  
  
  ### RMSE ###
  result.A3.rmse[i,1] <- sqrt(mean((predy-apply(x,1,output.A3high))^2)) # closed form
  result.A3.rmse[i,2] <- sqrt(mean((pred.muficokm$mean-apply(x,1,output.A3high))^2)) # Cokriging
  
  result.A3.meancrps[i,1] <- mean(crps(apply(x,1,output.A3high), predy, predsig2)) # closed form
  result.A3.meancrps[i,2] <- mean(crps(apply(x,1,output.A3high), pred.muficokm$mean, pred.muficokm$sig2)) # Cokriging

  result.A3.comptime[i,1] <- toc.closed - tic.closed
  result.A3.comptime[i,2] <- toc.cokm - tic.cokm
}


par(mfrow=c(1,1))
#RMSE comparison#
apply(result.A3.rmse, 2, mean) # 0.10, 0.13, 0.13, 0.05
table(apply(result.A3.rmse, 1, which.min))
boxplot(result.A3.rmse)

#CRPS comparison, The smaller, the better
apply(result.A3.meancrps, 2, mean)
table(apply(result.A3.meancrps, 1, which.min))
boxplot(result.A3.meancrps)




