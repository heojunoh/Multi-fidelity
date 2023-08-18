library(lhs)
library(laGP)

costmatc2 <- list(NA)
rmsematc2 <- list(NA)
### synthetic function ###
fl <- function(x, l){
  term1 <- sin(2*pi*x)
  term2 <- 0.2 * sin(8*pi*x)
  
  term1 + term2*5*0.8^l + (term1+term2)^3 + exp(-2*term1*term2)
}


### training data ###
n1 <- 12; n2 <- 8; n3 <- 5

for(kk in 1:10){
  set.seed(kk)
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
  
  ### test data ###
  x <- seq(0,1,0.01)
  
  
  ### closed ###
  fit.closed <- RNAmf2(X1, y1, X2, y2, X3, y3, kernel="sqex", constant=TRUE)
  predy <- predRNAmf2(fit.closed, x)$mu
  predsig2 <- predRNAmf2(fit.closed, x)$sig2
  
  ### RMSE ###
  sqrt(mean((predy-fl(x, l=5))^2)) # closed form
  
  
  ### IMSPE ###
  Icurrent <- mean(predsig2) # current IMSPE
  Icurrent
  
  ### Add 1 points and calculate IMSPE ###
  which.max(predsig2)
  
  intgvr <- integvar2(x[which.max(predsig2)], fit.closed, mc.sample=5)
  
  Icand1fast <- intgvr$intvar1
  Icand2fast <- intgvr$intvar2
  Icand3fast <- intgvr$intvar3
  
  ### Fast update; Equation 6.6. in Surrogates ###
  ### ALC; How much can be improved. Equation 6.6. in Surrogates ###
  alcfast <- c(Icurrent - Icand1fast, Icurrent - Icand2fast, Icurrent - Icand3fast)
  alcfast
  
  ### cost; 1, 2, 3 ###
  which.max(alcfast/c(2,(2+6),(2+6+12)))
  alcfast/c(2,(2+6),(2+6+12))
  
  
  chosen <- matrix(0, ncol=2)
  chosen[1,1] <- which.max(alcfast/c(2,(2+6),(2+6+12)))
  chosen[1,2] <- which.max(predsig2)
  
  
  ### Plotting the chosen point ###
  plot(x, predy, type="l", lwd=2, col=3, ylim=c(0,5))
  lines(x, predy+1.96*sqrt(predsig2*length(y2)/(length(y2)-2)), col=3, lty=2)
  lines(x, predy-1.96*sqrt(predsig2*length(y2)/(length(y2)-2)), col=3, lty=2)
  
  curve(fl(x, l=5),add=TRUE, col=1,lwd=2,lty=2) # high fidelity(TRUE); Black
  
  points(X1, y1, pch="1", col="red")
  points(X2, y2, pch="2", col="red")
  points(X3, y3, pch="3", col="red")
  
  text(x[which.min(Icand1fast)], predy[which.min(Icand1fast)], expression("1*"), col="red")
  text(x[which.min(Icand2fast)], predy[which.min(Icand2fast)], expression("2*"), col="red")
  text(x[which.min(Icand3fast)], predy[which.min(Icand3fast)], expression("3*"), col="red")
  
  nonlinear.cost <- 0
  nonlinear.error <- sqrt(mean((predy-fl(x, l=5))^2))
  
  Iselect <- IMSPEselect3(x[chosen[nrow(chosen),2]], fit.closed, level=chosen[nrow(chosen),1])
  
  
  
  #################
  ### Add point ###
  #################
  while(nonlinear.cost[length(nonlinear.cost)] < 100){ # if total cost is less than the budget
    
    ### closed ###
    predy <- predRNAmf2(Iselect$fit, x)$mu
    predsig2 <- predRNAmf2(Iselect$fit, x)$sig2
    
    ### RMSE ###
    nonlinear.error <- c(nonlinear.error, sqrt(mean((predy-fl(x, l=5))^2))) # closed form
    if(chosen[nrow(chosen),1] == 1){
      nonlinear.cost[length(nonlinear.cost)+1] <- nonlinear.cost[length(nonlinear.cost)]+2
    }else if(chosen[nrow(chosen),1] == 2){
      nonlinear.cost[length(nonlinear.cost)+1] <- nonlinear.cost[length(nonlinear.cost)]+(2+6)
    }else{
      nonlinear.cost[length(nonlinear.cost)+1] <- nonlinear.cost[length(nonlinear.cost)]+(2+6+12)
    }
    
    #############
    ### IMSPE ###
    #############
    Icurrent <- mean(predsig2)
    
    ### Add 1 points to the low-fidelity data ###
    intgvr <- integvar2(x[which.max(predsig2)], fit.closed, mc.sample=5)
    
    Icand1fast <- intgvr$intvar1
    Icand2fast <- intgvr$intvar2
    Icand3fast <- intgvr$intvar3
    
    ### Fast update; Equation 6.6. in Surrogates ###
    ### ALC; How much can be improved. Equation 6.6. in Surrogates ###
    alcfast <- c(Icurrent - Icand1fast, Icurrent - Icand2fast, Icurrent - Icand3fast)
    alcfast
    
    ### cost; 1, 2, 3 ###
    which.max(alcfast/c(2,(2+6),(2+6+12)))
    alcfast/c(2,(2+6),(2+6+12))
    
    
    chosen <- rbind(chosen, c(which.max(alcfast/c(2,(2+6),(2+6+12))), which.max(predsig2)))
    Iselect <- IMSPEselect3(x[chosen[nrow(chosen),2]], Iselect$fit, level=chosen[nrow(chosen),1])
    
    if(nonlinear.cost[length(nonlinear.cost)] >= 100){break}
    
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
}
costmatc2
rmsematc2



