library(lhs)
library(laGP)
library(plgp)
library(MuFiCokriging)

costmatc <- list(NA)
rmsematc <- list(NA)
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
  fit.closed <- RNAmf(X1, y1, X2, y2, kernel="sqex", constant=TRUE)
  predy <- predRNAmf(fit.closed, x)$mu
  predsig2 <- predRNAmf(fit.closed, x)$sig2

  ### RMSE ###
  sqrt(mean((predy-f2(x))^2)) # closed form


  ### IMSPE ###
  Icurrent <- mean(predsig2) # current IMSPE
  Icurrent

  ### Add 1 points and calculate IMSPE ###
  intgvr <- integvar(x, fit.closed, mc.sample=10)

  Icand1fast <- intgvr$intvar1
  Icand2fast <- intgvr$intvar2

  plot(x, Icand1fast, type="l", lwd=2, col=3, ylim=range(Icand1fast))
  plot(x, Icand2fast, type="l", lwd=2, col=3, ylim=range(Icand2fast))

  which.min(Icand1fast)
  which.min(Icand2fast)

  ### Fast update; Equation 6.6. in Surrogates ###
  ### ALC; How much can be improved. Equation 6.6. in Surrogates ###
  alcfast <- c(Icurrent - Icand1fast[which.min(Icand1fast)], Icurrent - Icand2fast[which.min(Icand2fast)])
  alcfast

  ### cost; 1, 2, 3 ###
  which.max(alcfast/c(1,(1+3)))
  alcfast/c(1,(1+3))


  chosen <- matrix(0, ncol=2)
  chosen[1,1] <- which.max(alcfast/c(1,(1+3)))
  chosen[1,2] <- which.min(cbind(Icand1fast, Icand2fast)[,chosen[1,1]])


  ### Plotting the chosen point ###
  plot(x, predy, type="l", lwd=2, col=3, ylim=c(-2,1))
  lines(x, predy+1.96*sqrt(predsig2*length(y2)/(length(y2)-2)), col=3, lty=2)
  lines(x, predy-1.96*sqrt(predsig2*length(y2)/(length(y2)-2)), col=3, lty=2)

  curve(f2(x),add=TRUE, col=1,lwd=2,lty=2) # high fidelity(TRUE); Black

  points(X1, y1, pch="1", col="red")
  points(X2, y2, pch="2", col="red")

  text(x[which.min(Icand1fast)], predy[which.min(Icand1fast)], expression("1*"), col="red")
  text(x[which.min(Icand2fast)], predy[which.min(Icand2fast)], expression("2*"), col="red")

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
      nonlinear.cost[length(nonlinear.cost)+1] <- nonlinear.cost[length(nonlinear.cost)]+1
    }else{
      nonlinear.cost[length(nonlinear.cost)+1] <- nonlinear.cost[length(nonlinear.cost)]+(1+3)
    }

    #############
    ### IMSPE ###
    #############
    Icurrent <- mean(predsig2)

    ### Add 1 points to the low-fidelity data ###
    intgvr <- integvar(x, Iselect$fit, mc.sample=10)

    Icand1fast <- intgvr$intvar1
    Icand2fast <- intgvr$intvar2

    which.min(Icand1fast)
    which.min(Icand2fast)

    ### Fast update; Equation 6.6. in Surrogates ###
    ### ALC; How much can be improved. Equation 6.6. in Surrogates ###
    alcfast <- c(Icurrent - Icand1fast[which.min(Icand1fast)], Icurrent - Icand2fast[which.min(Icand2fast)])
    alcfast

    ### cost; 1, 2, 3 ###
    which.max(alcfast/c(1,(1+3)))
    alcfast/c(1,(1+3))


    chosen <- rbind(chosen, c(which.max(alcfast/c(1,(1+3))), which.min(cbind(Icand1fast, Icand2fast)[,which.max(alcfast/c(1,(1+3)))])))
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

  text(x[which.min(Icand1fast)], predy[which.min(Icand1fast)], expression("1*"), col="red")
  text(x[which.min(Icand2fast)], predy[which.min(Icand2fast)], expression("2*"), col="red")


  ### Save results ###
  costmatc[[kk]] <- nonlinear.cost
  rmsematc[[kk]] <- nonlinear.error
}
costmatc
rmsematc



