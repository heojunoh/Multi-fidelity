library(lhs)
library(laGP)

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
n1 <- 12; n2 <- 9

for(kk in 1:10){
  set.seed(kk)
  X1 <- maximinLHS(n1, 1)
  X2 <- maximinLHS(n2, 1)
  
  NestDesign <- NestedDesignBuild(design = list(X1,X2))
  
  X1 <- NestDesign$PX
  X2 <- ExtractNestDesign(NestDesign,2)
  
  y1 <- f1(X1)
  y2 <- f2(X2)
  
  ### test data ###
  x <- seq(0,1,0.01)
  
  
  ### RNAmf ###
  fit.RNAmf <- RNAmf(X1, y1, X2, y2, kernel="sqex", constant=TRUE)
  predy <- predRNAmf(fit.RNAmf, x)$mu
  predsig2 <- predRNAmf(fit.RNAmf, x)$sig2
  
  ### RMSE ###
  sqrt(mean((predy-f2(x))^2)) # closed form
  
  
  ### current uncertainty at each level ###
  unc <- c(max(pred.GP(fit.RNAmf$fit1, x)$sig2), max(predsig2))
  uncind <- c(which.max(pred.GP(fit.RNAmf$fit1, x)$sig2), which.max(predsig2))
  
  ### cost; 1, 2, 3 ###
  which.max(unc/c(2,(2+8)))
  unc/c(2,(2+8))
  
  
  chosen <- matrix(0, ncol=2)
  chosen[1,1] <- which.max(unc/c(2,(2+8)))
  chosen[1,2] <- uncind[which.max(unc/c(2,(2+8)))]
  
  
  ### Plotting the chosen point ###
  plot(x, predy, type="l", lwd=2, col=3, ylim=c(-2,1))
  lines(x, predy+1.96*sqrt(predsig2*length(y2)/(length(y2)-2)), col=3, lty=2)
  lines(x, predy-1.96*sqrt(predsig2*length(y2)/(length(y2)-2)), col=3, lty=2)
  
  curve(f2(x),add=TRUE, col=1,lwd=2,lty=2) # high fidelity(TRUE); Black
  
  points(X1, y1, pch="1", col="red")
  points(X2, y2, pch="2", col="red")
  
  text(x[uncind[1]], predy[uncind[1]], expression("1*"), col="red")
  text(x[uncind[2]], predy[uncind[2]], expression("2*"), col="red")
  
  nonlinear.cost <- 0
  nonlinear.error <- sqrt(mean((predy-f2(x))^2))
  
  Iselect <- IMSPEselect1(x[chosen[nrow(chosen),2]], fit.RNAmf, level=chosen[nrow(chosen),1])
  
  
  #################
  ### Add point ###
  #################
  while(nonlinear.cost[length(nonlinear.cost)] < 100){ # if total cost is less than the budget
    
    ### RNAmf ###
    predy <- predRNAmf(Iselect$fit, x)$mu
    predsig2 <- predRNAmf(Iselect$fit, x)$sig2
    
    ### RMSE ###
    nonlinear.error <- c(nonlinear.error, sqrt(mean((predy-f2(x))^2))) # closed form
    if(chosen[nrow(chosen),1] == 1){
      nonlinear.cost[length(nonlinear.cost)+1] <- nonlinear.cost[length(nonlinear.cost)]+2
    }else{
      nonlinear.cost[length(nonlinear.cost)+1] <- nonlinear.cost[length(nonlinear.cost)]+(2+8)
    }
    
    #############
    ### IMSPE ###
    #############
    
    ### current uncertainty at each level ###
    unc <- c(max(pred.GP(Iselect$fit$fit1, x)$sig2), max(predsig2))
    uncind <- c(which.max(pred.GP(Iselect$fit$fit1, x)$sig2), which.max(predsig2))
    
    ### cost; 1, 2, 3 ###
    which.max(unc/c(2,(2+8)))
    unc/c(2,(2+8))
    
    
    chosen <- rbind(chosen, c(which.max(unc/c(2,(2+8))), uncind[which.max(unc/c(2,(2+8)))]))
    Iselect <- IMSPEselect1(x[chosen[nrow(chosen),2]], Iselect$fit, level=chosen[nrow(chosen),1])
    
    if(nonlinear.cost[length(nonlinear.cost)] >= 100){break}
    
  }
  
  
  ### Plotting the chosen point ###
  plot(x, predy, type="l", lwd=2, col=3,
       ylim=c(-2,1))
  lines(x, predy+1.96*sqrt(predsig2*length(y2)/(length(y2)-2)), col=3, lty=2)
  lines(x, predy-1.96*sqrt(predsig2*length(y2)/(length(y2)-2)), col=3, lty=2)
  
  curve(f2(x),add=TRUE, col=1,lwd=2,lty=2) # high fidelity(TRUE); Black
  
  points(t(t(Iselect$fit$fit1$X)*attr(Iselect$fit$fit1$X,"scaled:scale")+attr(Iselect$fit$fit1$X,"scaled:center")), Iselect$fit$fit1$y, pch="1", col="red")
  points(t(t(Iselect$fit$fit2$X)*attr(Iselect$fit$fit2$X,"scaled:scale")+attr(Iselect$fit$fit2$X,"scaled:center"))[,1], Iselect$fit$fit2$y, pch="2", col="red")
  
  text(x[uncind[1]], predy[uncind[1]], expression("1*"), col="red")
  text(x[uncind[2]], predy[uncind[2]], expression("2*"), col="red")
  
  
  ### Save results ###
  costmatc3[[kk]] <- nonlinear.cost
  rmsematc3[[kk]] <- nonlinear.error
}
costmatc3
rmsematc3



