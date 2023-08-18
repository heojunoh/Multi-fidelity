library(lhs)
library(laGP)
source("GP.R")
source("KOH.R")
source("closed.R")
source("IMSPE1.R")
source("IMSPE2.R")
source("IMSPE3.R")

fit.KOH <- function(X1, X2, Y1, Y2, g=eps){ # need to change function for another example
  
  ### KOH method ###
  Y1d2 <- f1(X2)
  
  ### estimating first order ###
  fit.KOHGP1 <- KOHGP(X1, Y1)
  b1 <- 1/fit.KOHGP1$theta
  sig2_1 <- fit.KOHGP1$tau2hat
  
  ### estimating second order ###
  # KOH(X2, Y2, Y1d2)
  rho1 <- KOH(X2, Y2, Y1d2)$rho
  b2 <- 1/KOH(X2, Y2, Y1d2)$theta
  sig2_2 <- KOH(X2, Y2, Y1d2)$tau2hat
  
  return(list(b=c(b1, b2), rho=rho1, tau2hat=c(sig2_1, sig2_2), g=g, X1=X1, X2=X2, Y1=Y1, Y2=Y2))
}

pred.KOH <- function(fit, x){ # need to change function for another example
  
  X1 <- fit$X1
  X2 <- fit$X2
  Y1 <- fit$Y1
  Y2 <- fit$Y2
  
  b <- fit$b
  rho <- fit$rho
  tau2hat <- fit$tau2hat
  g <- fit$g
  
  
  ### prediction of 2nd order KOH ###
  tx1 <- cbind(rho*tau2hat[1]*covar.sep(x, X1, d=1/b[1], g=g), 
               rho^2*tau2hat[1]*covar.sep(x, X2, d=1/b[1], g=g) + tau2hat[2]*covar.sep(x, X2, d=1/b[2], g=g))
  
  V1 <- tau2hat[1]*covar.sep(X1, d=1/b[1], g=g)
  V12 <- rho*tau2hat[1]*covar.sep(X1, X2, d=1/b[1], g=0)
  V2 <- rho^2*tau2hat[1]*covar.sep(X2, d=1/b[1], g=g) + tau2hat[2]*covar.sep(X2, d=1/b[2], g=g)
  
  V_2 <- rbind(cbind(V1, V12), cbind(t(V12), V2))
  
  mx1 <- tx1 %*% solve(V_2) %*% c(Y1, Y2)
  
  ### posterior variance ###
  koh.var1 <- pmax(0, diag(tau2hat[2]*covar.sep(as.matrix(x), d=1/b[2], g=g) + tau2hat[1]*rho^2*covar.sep(as.matrix(x), d=1/b[1], g=g) - tx1 %*% solve(V_2)%*%t(tx1)))
  
  return(list(mu=mx1, sig2=koh.var1))
}

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
  Icand1fast <- c(rep(0, length(g))) # IMSPE candidates
  Icand2fast <- c(rep(0, length(g))) # IMSPE candidates
  
  for(i in 1:length(Icand1fast)){ # no true, no need to fit just pred
    Icand1fast[i] <- IMSPEfast1(g, g[i], fit.closed, level=1)$IMSPE
  }
  for(i in 1:length(Icand2fast)){ # no true, no need to fit just pred
    Icand2fast[i] <- IMSPEfast1(g, g[i], fit.closed, level=2)$IMSPE
  }
  
  # plot(g, Icand1fast, type="l", lwd=2, col=3, ylim=range(Icand1fast))
  # plot(g, Icand2fast, type="l", lwd=2, col=3, ylim=range(Icand2fast))
  
  which.min(Icand1fast)
  which.min(Icand2fast)
  
  ### Fast update; Equation 6.6. in Surrogates ###
  ### ALC; How much can be improved. Equation 6.6. in Surrogates ###
  alcfast <- c(Icurrent - Icand1fast[which.min(Icand1fast)], Icurrent - Icand2fast[which.min(Icand2fast)])
  alcfast
  
  ### cost; 1, 2, 3 ###
  which.max(alcfast/c(2,(2+8)))
  alcfast/c(2,(2+8))
  
  
  chosen <- matrix(0, ncol=2)
  chosen[1,1] <- which.max(alcfast/c(2,(2+8)))
  chosen[1,2] <- which.min(cbind(Icand1fast, Icand2fast)[,chosen[1,1]])
  
  
  ### Plotting the chosen point ###
  plot(x, predy, type="l", lwd=2, col=3, ylim=c(-2,1)) 
  lines(x, predy+1.96*sqrt(predsig2*length(y2)/(length(y2)-2)), col=3, lty=2)
  lines(x, predy-1.96*sqrt(predsig2*length(y2)/(length(y2)-2)), col=3, lty=2)
  
  curve(f2(x),add=TRUE, col=1,lwd=2,lty=2) # high fidelity(TRUE); Black
  
  points(X1, y1, pch="1", col="red")
  points(X2, y2, pch="2", col="red")
  
  text(g[which.min(Icand1fast)], predy[which.min(Icand1fast)], expression("1*"), col="red")
  text(g[which.min(Icand2fast)], predy[which.min(Icand2fast)], expression("2*"), col="red")
  
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
    Icand1fast <- c(rep(0, length(g))) # IMSPE candidates
    Icand2fast <- c(rep(0, length(g))) # IMSPE candidates
    
    for(i in 1:length(Icand1fast)){ # no true, no need to fit just pred
      if(any(chosen[,2]==i)){Icand1fast[i] <- 0}else{
        Icand1fast[i] <- IMSPEfast1(g, g[i], Iselect$fit, level=1)$IMSPE
      }
    }
    for(i in 1:length(Icand2fast)){ # no true, no need to fit just pred
      if(any(chosen[,2]==i)){Icand2fast[i] <- 0}else{
        Icand2fast[i] <- IMSPEfast1(g, g[i], Iselect$fit, level=2)$IMSPE
      }
    }
    
    if(any(Icand1fast==0)){Icand1fast[which(Icand1fast==0)] <-  max(Icand1fast)}
    if(any(Icand2fast==0)){Icand2fast[which(Icand2fast==0)] <-  max(Icand2fast)}
    
    which.min(Icand1fast)
    which.min(Icand2fast)
    
    ### Fast update; Equation 6.6. in Surrogates ###
    ### ALC; How much can be improved. Equation 6.6. in Surrogates ###
    alcfast <- c(Icurrent - Icand1fast[which.min(Icand1fast)], Icurrent - Icand2fast[which.min(Icand2fast)])
    alcfast
    
    ### cost; 1, 2, 3 ###
    which.max(alcfast/c(2,(2+8)))
    alcfast/c(2,(2+8))
    
    
    chosen <- rbind(chosen, c(which.max(alcfast/c(2,(2+8))), which.min(cbind(Icand1fast, Icand2fast)[,which.max(alcfast/c(2,(2+8)))])))
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
  
  text(g[which.min(Icand1fast)], predy[which.min(Icand1fast)], expression("1*"), col="red")
  text(g[which.min(Icand2fast)], predy[which.min(Icand2fast)], expression("2*"), col="red")
  
  
  ### Save results ###
  costmatc[[kk]] <- nonlinear.cost
  rmsematc[[kk]] <- nonlinear.error
}
costmatc
rmsematc



