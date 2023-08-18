library(lhs)
library(laGP)
source("GP.R")
source("KOH.R")
source("closed.R")

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

costmatk <- list(NA)
rmsematk <- list(NA)
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
  
  ### KOH method ###
  fit.KOH2 <- fit.KOH(X1, X2, y1, y2)
  pred.KOH2 <- pred.KOH(fit.KOH2, x)
  mx1 <- pred.KOH2$mu
  koh.var1 <- pred.KOH2$sig2
  
  ### RMSE ###
  sqrt(mean((mx1-f2(x))^2)) # KOH
  
  
  ### IMSPE ###
  g <- seq(0,1,0.01) # more than 1-dim, use expand.grid()
  Icurrent <- mean(koh.var1) # current IMSPE
  Icurrent
  
  ### Add 1 points and calculate IMSPE ###
  IcandKOH1 <- c(rep(0, length(g))) # IMSPE candidates
  IcandKOH2 <- c(rep(0, length(g))) # IMSPE candidates
  
  for(i in 1:length(IcandKOH1)){ # no true, no need to fit just pred
    IcandKOH1[i] <- IMSPEKOH1(g, g[i], fit.KOH2, level=1)
  }
  for(i in 1:length(IcandKOH2)){ # no true, no need to fit just pred
    IcandKOH2[i] <- IMSPEKOH1(g, g[i], fit.KOH2, level=2)
  }
  
  # plot(g, IcandKOH1, type="l", lwd=2, col=3, ylim=range(IcandKOH1))
  # plot(g, IcandKOH2, type="l", lwd=2, col=3, ylim=range(IcandKOH2))
  
  which.min(IcandKOH1)
  which.min(IcandKOH2)
  
  ### Fast update; Equation 6.6. in Surrogates ###
  ### ALC; How much can be improved. Equation 6.6. in Surrogates ###
  alcfast <- c(Icurrent - IcandKOH1[which.min(IcandKOH1)], Icurrent - IcandKOH2[which.min(IcandKOH2)])
  alcfast
  
  ### cost; 1, 2, 3 ###
  which.max(alcfast/c(2,8))
  alcfast/c(2,8)
  
  
  chosen <- matrix(0, ncol=2)
  chosen[1,1] <- which.max(alcfast/c(2,8))
  chosen[1,2] <- which.min(cbind(IcandKOH1, IcandKOH2)[,chosen[1,1]])
  
  
  ### Plotting the chosen point ###
  plot(x, mx1, type="l", lwd=2, col=3, ylim=c(-2,1))
  lines(x, mx1+1.96*sqrt(koh.var1*length(y2)/(length(y2)-2)), col=3, lty=2)
  lines(x, mx1-1.96*sqrt(koh.var1*length(y2)/(length(y2)-2)), col=3, lty=2)
  
  curve(f2(x),add=TRUE, col=1,lwd=2,lty=2) # high fidelity(TRUE); Black
  
  points(X1, y1, pch="1", col="red")
  points(X2, y2, pch="2", col="red")
  
  text(g[which.min(IcandKOH1)], mx1[which.min(IcandKOH1)], expression("1*"), col="red")
  text(g[which.min(IcandKOH2)], mx1[which.min(IcandKOH2)], expression("2*"), col="red")
  
  nonlinear.cost <- 0
  nonlinear.error <- sqrt(mean((mx1-f2(x))^2))
  
  Iselect <- IMSPEKOHselect1(g, g[chosen[nrow(chosen),2]], fit.KOH2, level=chosen[nrow(chosen),1]) 
  
  
  #################
  ### Add point ###
  #################
  while(nonlinear.cost[length(nonlinear.cost)] < 100){ # if total cost is less than the budget
    
    ### predictive ###
    pred.KOH2 <- pred.KOH(Iselect, x)
    mx1 <- pred.KOH2$mu
    koh.var1 <- pred.KOH2$sig2
    
    ### RMSE ###  
    nonlinear.error <- c(nonlinear.error, sqrt(mean((mx1-f2(x))^2))) # closed form
    if(chosen[nrow(chosen),1] == 1){
      nonlinear.cost[length(nonlinear.cost)+1] <- nonlinear.cost[length(nonlinear.cost)]+2
    }else{
      nonlinear.cost[length(nonlinear.cost)+1] <- nonlinear.cost[length(nonlinear.cost)]+8
    }
    
    #############
    ### IMSPE ###
    #############
    g <- seq(0,1,0.01) # more than 1-dim, use expand.grid()
    Icurrent <- mean(koh.var1)
    
    ### Add 1 points and calculate IMSPE ###
    IcandKOH1 <- c(rep(0, length(g))) # IMSPE candidates
    IcandKOH2 <- c(rep(0, length(g))) # IMSPE candidates
    
    for(i in 1:length(IcandKOH1)){ # no true, no need to fit just pred
      if(any(chosen[,2]==i)){IcandKOH1[i] <- 0}else{
        IcandKOH1[i] <- IMSPEKOH1(g, g[i], Iselect, level=1)
      }
    }
    for(i in 1:length(IcandKOH2)){ # no true, no need to fit just pred
      if(any(chosen[,2]==i)){IcandKOH2[i] <- 0}else{
        IcandKOH2[i] <- IMSPEKOH1(g, g[i], Iselect, level=2)
      }
    }
    
    if(any(IcandKOH1==0)){IcandKOH1[which(IcandKOH1==0)] <-  max(IcandKOH1)}
    if(any(IcandKOH2==0)){IcandKOH2[which(IcandKOH2==0)] <-  max(IcandKOH2)}
    
    which.min(IcandKOH1)
    which.min(IcandKOH2)
    
    ### Fast update; Equation 6.6. in Surrogates ###
    ### ALC; How much can be improved. Equation 6.6. in Surrogates ###
    alcfast <- c(Icurrent - IcandKOH1[which.min(IcandKOH1)], Icurrent - IcandKOH2[which.min(IcandKOH2)])
    alcfast
    
    ### cost; 1, 2, 3 ###
    which.max(alcfast/c(2,8))
    alcfast/c(2,8)
    
    
    chosen <- rbind(chosen, c(which.max(alcfast/c(2,8)), which.min(cbind(IcandKOH1, IcandKOH2)[,which.max(alcfast/c(2,8))])))
    Iselect <- IMSPEKOHselect1(g, g[chosen[nrow(chosen),2]], Iselect, level=chosen[nrow(chosen),1])
    
    if(nonlinear.cost[length(nonlinear.cost)] >= 100){break}
    
  }
  
  
  ### Plotting the chosen point ###
  plot(x, mx1, type="l", lwd=2, col=3,
       ylim=c(-2,1))
  lines(x, mx1+1.96*sqrt(koh.var1*length(y2)/(length(y2)-2)), col=3, lty=2)
  lines(x, mx1-1.96*sqrt(koh.var1*length(y2)/(length(y2)-2)), col=3, lty=2)
  
  curve(f2(x),add=TRUE, col=1,lwd=2,lty=2) # high fidelity(TRUE); Black
  
  points(Iselect$X1, Iselect$Y1, pch="1", col="red") 
  points(Iselect$X2, Iselect$Y2, pch="2", col="red")
  
  text(g[which.min(IcandKOH1)], mx1[which.min(IcandKOH1)], expression("1*"), col="red")
  text(g[which.min(IcandKOH2)], mx1[which.min(IcandKOH2)], expression("2*"), col="red")
  
  
  ### Save results ###
  costmatk[[kk]] <- nonlinear.cost
  rmsematk[[kk]] <- nonlinear.error
}
costmatk
rmsematk



