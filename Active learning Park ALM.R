library(lhs)
library(laGP)
source("GP.R")
source("closed.R")
source("IMSPE1.R")
source("IMSPE2.R")
source("IMSPE3.R")


costmatc2 <- list(NA)
rmsematc2 <- list(NA)
### synthetic function ###
park91a <- function(xx)
{
  x1 <- xx[1]
  x2 <- xx[2]
  x3 <- xx[3]
  x4 <- xx[4]
  
  term1a <- x1 / 2
  term1b <- sqrt(1 + (x2+x3^2)*x4/(x1^2)) - 1
  term1 <- term1a * term1b
  
  term2a <- x1 + 3*x4
  term2b <- exp(1 + sin(x3))
  term2 <- term2a * term2b
  
  y <- term1 + term2
  return(y)
}

park91alc <- function(xx)
{
  x1 <- xx[1]
  x2 <- xx[2]
  x3 <- xx[3]
  x4 <- xx[4]
  
  yh <- park91a(xx)
  
  term1 <- (1+sin(x1)/10) * yh
  term2 <- -2*x1 + x2^2 + x3^2
  
  y <- term1 + term2 + 0.5
  return(y)
}

### training data ###
n1 <- 15; n2 <- 8
d <- 4
for(kk in 1:10){
  set.seed(kk)
  X1 <- maximinLHS(n1, d)
  X2 <- maximinLHS(n2, d)
  
  NestDesign <- NestedDesignBuild(design = list(X1,X2))
  
  X1 <- NestDesign$PX
  X2 <- ExtractNestDesign(NestDesign,2)
  
  y1 <- apply(X1,1,park91alc)
  y2 <- apply(X2,1,park91a)
  
  ### model fitting for f1 ###
  eps <- sqrt(.Machine$double.eps)
  fit.GP1 <- GP(X1, y1)
  
  ### model fitting using (x2, f1(x2)) ###
  w1.x2 <- pred.GP(fit.GP1, X2)$mu # can interpolate; nested
  X2new <- cbind(X2, w1.x2) # combine (X2, f1(x2))
  fit.GP2new <- GP(X2new, y2) # model fitting for f_M(X2, f1(x2))
  
  
  ### test data ###
  x <- maximinLHS(300, d)
  
  ### closed ###
  fit.closed <- closed(X1, y1, X2, y2, constant=TRUE)
  predy <- predclosed(fit.closed, x)$mu
  predsig2 <- predclosed(fit.closed, x)$sig2
  
  ### compared to single fidelity ###
  fit.GP2 <- GP(X2, y2)
  pred2 <- pred.GP(fit.GP2, x)
  
  ### direct fitting; not using closed form. f1(u) from (u, f1(u)) is random variable.
  x1.mu <- rnorm(nrow(x), mean=pred.GP(fit.GP1, x)$mu, sd=sqrt(pred.GP(fit.GP1, x)$sig2)) 
  xnew <- cbind(x, x1.mu) # Use mu of the input in the closed form
  pred2new <- pred.GP(fit.GP2new, xnew) # not closed form
  
  
  ### RMSE ###
  sqrt(sum((predy-apply(x,1,park91a))^2))/(sqrt(sum((apply(x,1,park91a))^2))) # closed form
  sqrt(sum((pred2new$mu-apply(x,1,park91a))^2))/(sqrt(sum((apply(x,1,park91a))^2))) # not closed form
  sqrt(sum((pred2$mu-apply(x,1,park91a))^2))/(sqrt(sum((apply(x,1,park91a))^2))) # single fidelity
  
  
  ### IMSPE ###
  Icurrent <- mean(predsig2) # current IMSPE
  Icurrent
  
  ### Add 1 points and calculate IMSPE ###
  Icand1fast <- IMSPEfast1(x, x[which.max(predsig2),], fit.closed, level=1)$IMSPE
  Icand2fast <- IMSPEfast1(x, x[which.max(predsig2),], fit.closed, level=2)$IMSPE
  
  # plot(g, Icand1fast, type="l", lwd=2, col=3, ylim=range(Icand1fast))
  # plot(g, Icand2fast, type="l", lwd=2, col=3, ylim=range(Icand2fast))
  
  which.min(Icand1fast)
  which.min(Icand2fast)
  
  ### Fast update; Equation 6.6. in Surrogates ###
  ### ALC; How much can be improved. Equation 6.6. in Surrogates ###
  alcfast <- c(Icurrent - Icand1fast, Icurrent - Icand2fast)
  alcfast
  
  ### cost; 1, 2, 3 ###
  which.max(alcfast/c(1,(1+10)))
  alcfast/c(1,(1+10))
  
  
  chosen <- matrix(0, ncol=2)
  chosen[1,1] <- which.max(alcfast/c(1,(1+10)))
  chosen[1,2] <- which.max(predsig2)
  
  
  nonlinear.cost <- 0
  nonlinear.error <- sqrt(sum((predy-apply(x,1,park91a))^2))/(sqrt(sum((apply(x,1,park91a))^2)))
  
  Iselect <- IMSPEselect1(x, x[chosen[nrow(chosen),2],], fit.closed, level=chosen[nrow(chosen),1]) 
  
  
  
  #################
  ### Add point ###
  #################
  while(nonlinear.cost[length(nonlinear.cost)] < 100){ # if total cost is less than the budget
    
    ### closed ###
    predy <- predclosed(Iselect$fit, x)$mu
    predsig2 <- predclosed(Iselect$fit, x)$sig2
    
    ### RMSE ###  
    nonlinear.error <- c(nonlinear.error, sqrt(sum((predy-apply(x,1,park91a))^2))/(sqrt(sum((apply(x,1,park91a))^2)))) # closed form
    if(chosen[nrow(chosen),1] == 1){
      nonlinear.cost[length(nonlinear.cost)+1] <- nonlinear.cost[length(nonlinear.cost)]+1
    }else{
      nonlinear.cost[length(nonlinear.cost)+1] <- nonlinear.cost[length(nonlinear.cost)]+(1+10)
    }
    
    #############
    ### IMSPE ###
    #############
    Icurrent <- Iselect$IMSPE
    
    ### Add 1 points to the low-fidelity data ###
    Icand1fast <- IMSPEfast1(x, x[which.max(predsig2),], Iselect$fit, level=1)$IMSPE
    Icand2fast <- IMSPEfast1(x, x[which.max(predsig2),], Iselect$fit, level=2)$IMSPE
    
    which.min(Icand1fast)
    which.min(Icand2fast)
    
    ### Fast update; Equation 6.6. in Surrogates ###
    ### ALC; How much can be improved. Equation 6.6. in Surrogates ###
    alcfast <- c(Icurrent - Icand1fast, Icurrent - Icand2fast)
    alcfast
    
    ### cost; 1, 2, 3 ###
    which.max(alcfast/c(1,(1+10)))
    alcfast/c(1,(1+10))
    
    
    chosen <- rbind(chosen, c(which.max(alcfast/c(1,(1+10))), which.max(predsig2)))
    Iselect <- IMSPEselect1(x, x[chosen[nrow(chosen),2],], Iselect$fit, level=chosen[nrow(chosen),1])
    
    if(nonlinear.cost[length(nonlinear.cost)] >= 100){break}
    
  }
  
  
  ### Save results ###
  costmatc2[[kk]] <- nonlinear.cost
  rmsematc2[[kk]] <- nonlinear.error
}

costmatc2
rmsematc2



plot(costmatc[[1]], rmsematc[[1]], ylim=c(0, 0.1), lwd=2, type="l", col="green")
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
