library(lhs)
library(laGP)
library(plgp)
library(MuFiCokriging)

costmatc3 <- list(NA)
rmsematc3 <- list(NA)
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
  
  Iselect <- ALMC_two_level(x, fit.closed, 100, c(1,3), list(f1, f2))
  
  
  #################
  ### Add point ###
  #################
  while(nonlinear.cost[length(nonlinear.cost)] < 100){ # if total cost is less than the budget
    ### closed ###
    predy <- predRNAmf(Iselect$fit, x)$mu
    predsig2 <- predRNAmf(Iselect$fit, x)$sig2
    
    ### RMSE ###
    nonlinear.error <- c(nonlinear.error, sqrt(mean((predy-f2(x))^2))) # closed form
    if(Iselect$chosen$level == 1){
      nonlinear.cost[length(nonlinear.cost)+1] <- nonlinear.cost[length(nonlinear.cost)]+1
    }else{
      nonlinear.cost[length(nonlinear.cost)+1] <- nonlinear.cost[length(nonlinear.cost)]+(1+3)
    }
    print(nonlinear.cost[length(nonlinear.cost)])
    print(nonlinear.error[length(nonlinear.error)])
    
    if(nonlinear.cost[length(nonlinear.cost)] >= 100){break}

    ### update the next point ###
    Iselect <- ALMC_two_level(x, Iselect$fit, 100, c(1,3), list(f1, f2))
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
  # text(x[uncind[1]], predy[uncind[1]], expression("1*"), col="red")
  # text(x[uncind[2]], predy[uncind[2]], expression("2*"), col="red")
  
  
  ### Save results ###
  costmatc3[[kk]] <- nonlinear.cost
  rmsematc3[[kk]] <- nonlinear.error
}
costmatc3
rmsematc3



