library(lhs)
library(laGP)
source("GP.R")
source("KOH.R")
source("closed.R")
source("IMSPE1.R")
source("IMSPE2.R")
source("IMSPE3.R")

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
  
  ### model fitting for f1 ###
  eps <- sqrt(.Machine$double.eps)
  fit.GP1 <- GP(X1, y1, constant=TRUE)
  
  ### model fitting using (x2, f1(x2)) ###
  w1.x2 <- pred.GP(fit.GP1, X2)$mu # can interpolate; nested
  X2new <- cbind(X2, w1.x2) # combine (X2, f1(x2))
  fit.GP2new <- GP(X2new, y2, constant=TRUE) # model fitting for f_M(X2, f1(x2))
  
  ### test data ###
  x <- seq(0,1,0.01)
  
  
  ### closed ###
  fit.closed <- closed(X1, y1, X2, y2, constant=TRUE)
  predy <- predclosed(fit.closed, x)$mu
  predsig2 <- predclosed(fit.closed, x)$sig2
  
  
  ### compared to single fidelity ###
  fit.GP2 <- GP(X2, y2, constant=TRUE)
  pred2 <- pred.GP(fit.GP2, x)
  
  ### direct fitting; not using closed form. f1(u) and f_M(u) from (u, f_M(u, f1(u))) are random variables.
  w1.x <- c(rep(NA, length(x)))
  for(i in 1:length(x)){
    w1.x[i] <- mean(rnorm(10000, mean=pred.GP(fit.GP1, x[i])$mu, sd=sqrt(pred.GP(fit.GP1, x[i])$sig2)))
  }
  
  xxnew <- cbind(x, w1.x)
  pred2new <- pred.GP(fit.GP2new, xxnew) # not closed form
  
  
  ### RMSE ###
  sqrt(sum((predy-f2(x))^2))/(sqrt(sum((f2(x))^2))) # closed form
  sqrt(sum((pred2new$mu-f2(x))^2))/(sqrt(sum((f2(x))^2))) # not closed form
  sqrt(sum((pred2$mu-f2(x))^2))/(sqrt(sum((f2(x))^2))) # single fidelity
  
  
  ### IMSPE ###
  g <- seq(0,1,0.01) # more than 1-dim, use expand.grid()
  Icurrent <- mean(predsig2) # current IMSPE
  Icurrent
  
  ### Add 1 points and calculate IMSPE ###
  which.max(predsig2)
  
  Icand1fast <- IMSPEfast1(g, g[which.max(predsig2)], fit.closed, level=1)$IMSPE
  Icand2fast <- IMSPEfast1(g, g[which.max(predsig2)], fit.closed, level=2)$IMSPE
  
  # plot(g, Icand1fast, type="l", lwd=2, col=3, ylim=range(Icand1fast))
  # plot(g, Icand2fast, type="l", lwd=2, col=3, ylim=range(Icand2fast))
  
  Icand1fast
  Icand2fast
  
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
  
  text(g[which.max(predsig2)], predy[which.max(predsig2)], expression("1*"), col="red")
  text(g[which.max(predsig2)], predy[which.max(predsig2)], expression("2*"), col="red")
  
  nonlinear.cost <- 0
  nonlinear.error <- sqrt(sum((predy-f2(x))^2))/(sqrt(sum((f2(x))^2)))
  
  Iselect <- IMSPEselect1(g, g[chosen[nrow(chosen),2]], fit.closed, level=chosen[nrow(chosen),1]) 
  
  
  
  #################
  ### Add point ###
  #################
  while(nonlinear.cost[length(nonlinear.cost)] < 100){ # if total cost is less than the budget
    
    ### closed ###
    predy <- predclosed(Iselect$fit, x)$mu
    predsig2 <- predclosed(Iselect$fit, x)$sig2
    
    ### RMSE ###  
    nonlinear.error <- c(nonlinear.error, sqrt(sum((predy-f2(x))^2))/(sqrt(sum((f2(x))^2)))) # closed form
    if(chosen[nrow(chosen),1] == 1){
      nonlinear.cost[length(nonlinear.cost)+1] <- nonlinear.cost[length(nonlinear.cost)]+2
    }else{
      nonlinear.cost[length(nonlinear.cost)+1] <- nonlinear.cost[length(nonlinear.cost)]+(2+8)
    }
    
    #############
    ### IMSPE ###
    #############
    g <- seq(0,1,0.01) # more than 1-dim, use expand.grid()
    Icurrent <- Iselect$IMSPE
    
    ### Add 1 points to the low-fidelity data ###
    Icand1fast <- IMSPEfast1(g, g[which.max(predsig2)], fit.closed, level=1)$IMSPE
    Icand2fast <- IMSPEfast1(g, g[which.max(predsig2)], fit.closed, level=2)$IMSPE
    
    Icand1fast
    Icand2fast
    
    ### Fast update; Equation 6.6. in Surrogates ###
    ### ALC; How much can be improved. Equation 6.6. in Surrogates ###
    alcfast <- c(Icurrent - Icand1fast, Icurrent - Icand2fast)
    alcfast
    
    ### cost; 1, 2, 3 ###
    which.max(alcfast/c(2,(2+8)))
    alcfast/c(2,(2+8))
    
    
    chosen <- rbind(chosen, c(which.max(alcfast/c(2,(2+8))), which.max(predsig2)))
    Iselect <- IMSPEselect1(g, g[chosen[nrow(chosen),2]], Iselect$fit, level=chosen[nrow(chosen),1])
    
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
  
  text(g[which.max(predsig2)], predy[which.max(predsig2)], expression("1*"), col="red")
  text(g[which.max(predsig2)], predy[which.max(predsig2)], expression("2*"), col="red")
  
  
  ### Save results ###
  costmatc2[[kk]] <- nonlinear.cost
  rmsematc2[[kk]] <- nonlinear.error
}

costmatc2
rmsematc2



plot(costmatc[[1]], rmsematc[[1]], ylim=c(0, 0.65), lwd=2, type="l", col="green")
lines(costmatc[[2]], rmsematc[[2]], lwd=2, type="l", col="green")
lines(costmatc[[3]], rmsematc[[3]], lwd=2, type="l", col="green")
lines(costmatc[[4]], rmsematc[[4]], lwd=2, type="l", col="green")
lines(costmatc[[5]], rmsematc[[5]], lwd=2, type="l", col="green")
lines(costmatc[[6]], rmsematc[[6]], lwd=2, type="l", col="green")
lines(costmatc[[7]], rmsematc[[7]], lwd=2, type="l", col="green")
lines(costmatc[[8]], rmsematc[[8]], lwd=2, type="l", col="green")
lines(costmatc[[9]], rmsematc[[9]], lwd=2, type="l", col="green")
lines(costmatc[[10]], rmsematc[[10]], lwd=2, type="l", col="green")


lines(costmatc2[[1]], rmsematc2[[1]], lwd=2, type="l", col="red")
lines(costmatc2[[2]], rmsematc2[[2]], lwd=2, type="l", col="red")
lines(costmatc2[[3]], rmsematc2[[3]], lwd=2, type="l", col="red")
lines(costmatc2[[4]], rmsematc2[[4]], lwd=2, type="l", col="red")
lines(costmatc2[[5]], rmsematc2[[5]], lwd=2, type="l", col="red")
lines(costmatc2[[6]], rmsematc2[[6]], lwd=2, type="l", col="red")
lines(costmatc2[[7]], rmsematc2[[7]], lwd=2, type="l", col="red")
lines(costmatc2[[8]], rmsematc2[[8]], lwd=2, type="l", col="red")
lines(costmatc2[[9]], rmsematc2[[9]], lwd=2, type="l", col="red")
lines(costmatc2[[10]], rmsematc2[[10]], lwd=2, type="l", col="red")

legend("topright", c("ALC", "ALM"), lty=1, lwd=2, col=c("green", "red"))
