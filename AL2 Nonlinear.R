library(lhs)
library(laGP)
library(plgp)
library(MuFiCokriging)
library(doParallel)
library(foreach)
library(RNAmf)

eps <- sqrt(.Machine$double.eps)
crps <- function(x, mu, sig2){ # The smaller, the better (0 to infinity)
  if(any(sig2==0)) sig2[sig2==0] <- eps
  -sqrt(sig2)*(1/sqrt(pi)-2*dnorm((x-mu)/sqrt(sig2))-(x-mu)/sqrt(sig2)*(2*pnorm((x-mu)/sqrt(sig2))-1))
}

costmatc2 <- list(NA)
rmsematc2 <- list(NA)
crpsmatc2 <- list(NA)
time.each2 <- rep(0,10)
### synthetic function ###
f1 <- function(x)
{
  sin(8*pi*x)
}

f2 <- function(x)
{
  (x-sqrt(2))*(sin(8*pi*x))^2
}

### training data ###
n1 <- 13; n2 <- 8

for(kk in 1:10){
  time.start <- proc.time()[3]
  set.seed(kk)
  print(kk)
  X1 <- maximinLHS(n1, 1)
  X2 <- maximinLHS(n2, 1)

  NestDesign <- NestedDesignBuild(design = list(X1,X2))

  X1 <- NestDesign$PX
  X2 <- ExtractNestDesign(NestDesign,2)

  y1 <- f1(X1)
  y2 <- f2(X2)

  ### test data ###
  x <- seq(0,1,length.out=1000)


  ### closed ###
  fit.closed <- RNAmf(X1, y1, X2, y2, kernel="sqex", constant=TRUE)
  predy <- predRNAmf(fit.closed, x)$mu
  predsig2 <- predRNAmf(fit.closed, x)$sig2

  ### RMSE ###
  sqrt(mean((predy-f2(x))^2)) # closed form
  
  nonlinear.cost <- 0
  nonlinear.error <- sqrt(mean((predy-f2(x))^2))
  nonlinear.crps <- mean(crps(f2(x), predy, predsig2))
  
  Iselect <- ALC_two_level(x, fit.closed, 100, c(1,9), list(f1, f2), parallel=TRUE, ncore=10)


  #################
  ### Add point ###
  #################
  while(nonlinear.cost[length(nonlinear.cost)] < 100){ # if total cost is less than the budget
    ### closed ###
    predy <- predRNAmf(Iselect$fit, x)$mu
    predsig2 <- predRNAmf(Iselect$fit, x)$sig2

    ### RMSE ###
    nonlinear.error <- c(nonlinear.error, sqrt(mean((predy-f2(x))^2))) # RMSE
    nonlinear.crps <- c(nonlinear.crps, mean(crps(f2(x), predy, predsig2))) # CRPS
    if(Iselect$chosen$level == 1){
      nonlinear.cost[length(nonlinear.cost)+1] <- nonlinear.cost[length(nonlinear.cost)]+1
    }else{
      nonlinear.cost[length(nonlinear.cost)+1] <- nonlinear.cost[length(nonlinear.cost)]+(1+9)
    }
    print(nonlinear.cost[length(nonlinear.cost)])
    print(nonlinear.error[length(nonlinear.error)])
    
    if(nonlinear.cost[length(nonlinear.cost)] >= 100){break}
    
    ### update the next point ###
    Iselect <- ALC_two_level(x, Iselect$fit, 100, c(1,9), list(f1, f2), parallel=TRUE, ncore=10)
    save.image("/Users/junoh/Downloads/Perd AL2 1,9.RData")
  }


  # ### Plotting the chosen point ###
  # plot(x, predy, type="l", lwd=2, col=3,
  #      ylim=c(-2,1))
  # lines(x, predy+1.96*sqrt(predsig2*length(y2)/(length(y2)-2)), col=3, lty=2)
  # lines(x, predy-1.96*sqrt(predsig2*length(y2)/(length(y2)-2)), col=3, lty=2)
  # 
  # curve(f2(x),add=TRUE, col=1,lwd=2,lty=2) # high fidelity(TRUE); Black
  # 
  # points(t(t(Iselect$fit$fit1$X)*attr(Iselect$fit$fit1$X,"scaled:scale")+attr(Iselect$fit$fit1$X,"scaled:center")), Iselect$fit$fit1$y, pch="1", col="red")
  # points(t(t(Iselect$fit$fit2$X)*attr(Iselect$fit$fit2$X,"scaled:scale")+attr(Iselect$fit$fit2$X,"scaled:center"))[,1], Iselect$fit$fit2$y, pch="2", col="red")
  # 
  # text(x[which.max(predsig2)], predy[which.max(predsig2)], expression("1*"), col="red")
  # text(x[which.max(predsig2)], predy[which.max(predsig2)], expression("2*"), col="red")


  ### Save results ###
  costmatc2[[kk]] <- nonlinear.cost
  rmsematc2[[kk]] <- nonlinear.error
  crpsmatc2[[kk]] <- nonlinear.crps
  save.image("/Users/junoh/Downloads/Perd AL2 1,9.RData")
  
  time.each2[kk] <- proc.time()[3]- time.start
}
costmatc2
rmsematc2
crpsmatc2
time.each2



