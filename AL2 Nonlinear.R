library(lhs)
library(laGP)

costmatc2 <- list(NA)
rmsematc2 <- list(NA)
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


  ### closed ###
  fit.closed <- closed(X1, y1, X2, y2, kernel="sqex", constant=TRUE)
  predy <- predclosed(fit.closed, x)$mu
  predsig2 <- predclosed(fit.closed, x)$sig2

  ### RMSE ###
  sqrt(mean((predy-f2(x))^2)) # closed form


  ### IMSPE ###
  Icurrent <- mean(predsig2) # current IMSPE
  Icurrent

  ### Add 1 points and calculate IMSPE ###
  which.max(predsig2)

  intgvr <- integvar(x[which.max(predsig2)], fit.closed, mc.sample=100)

  Icand1fast <- intgvr$intvar1
  Icand2fast <- intgvr$intvar2

  ### Fast update; Equation 6.6. in Surrogates ###
  ### ALC; How much can be improved. Equation 6.6. in Surrogates ###
  alcfast <- c(Icurrent - Icand1fast, Icurrent - Icand2fast)
  alcfast

  ### cost; 1, 2, 3 ###
  which.max(alcfast/c(2,(2+8)))
  alcfast/c(2,(2+8))


  chosen <- matrix(0, ncol=2)
  chosen[1,1] <- which.max(alcfast/c(2,(2+8)))
  chosen[1,2] <- which.max(predsig2)


  ### Plotting the chosen point ###
  plot(x, predy, type="l", lwd=2, col=3, ylim=c(-2,1))
  lines(x, predy+1.96*sqrt(predsig2*length(y2)/(length(y2)-2)), col=3, lty=2)
  lines(x, predy-1.96*sqrt(predsig2*length(y2)/(length(y2)-2)), col=3, lty=2)

  curve(f2(x),add=TRUE, col=1,lwd=2,lty=2) # high fidelity(TRUE); Black

  points(X1, y1, pch="1", col="red")
  points(X2, y2, pch="2", col="red")

  text(x[which.max(predsig2)], predy[which.max(predsig2)], expression("1*"), col="red")
  text(x[which.max(predsig2)], predy[which.max(predsig2)], expression("2*"), col="red")

  nonlinear.cost <- 0
  nonlinear.error <- sqrt(mean((predy-f2(x))^2))

  Iselect <- IMSPEselect1(x[chosen[nrow(chosen),2]], fit.closed, level=chosen[nrow(chosen),1])



  #################
  ### Add point ###
  #################
  while(nonlinear.cost[length(nonlinear.cost)] < 100){ # if total cost is less than the budget

    ### closed ###
    predy <- predclosed(Iselect$fit, x)$mu
    predsig2 <- predclosed(Iselect$fit, x)$sig2

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
    Icurrent <- mean(predsig2)

    ### Add 1 points to the low-fidelity data ###
    intgvr <- integvar(x[which.max(predsig2)], fit.closed, mc.sample=100)

    Icand1fast <- intgvr$intvar1
    Icand2fast <- intgvr$intvar2

    ### Fast update; Equation 6.6. in Surrogates ###
    ### ALC; How much can be improved. Equation 6.6. in Surrogates ###
    alcfast <- c(Icurrent - Icand1fast, Icurrent - Icand2fast)
    alcfast

    ### cost; 1, 2, 3 ###
    which.max(alcfast/c(2,(2+8)))
    alcfast/c(2,(2+8))


    chosen <- rbind(chosen, c(which.max(alcfast/c(2,(2+8))), which.max(predsig2)))
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

  text(x[which.max(predsig2)], predy[which.max(predsig2)], expression("1*"), col="red")
  text(x[which.max(predsig2)], predy[which.max(predsig2)], expression("2*"), col="red")


  ### Save results ###
  costmatc2[[kk]] <- nonlinear.cost
  rmsematc2[[kk]] <- nonlinear.error
}
costmatc2
rmsematc2



