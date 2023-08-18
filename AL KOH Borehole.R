library(lhs)
library(laGP)
source("GP.R")
source("KOH.R")
source("closed.R")

fit.KOH <- function(X1, X2, Y1, Y2, g=eps){ # need to change function for another example
  
  ### KOH method ###
  Y1d2 <- apply(X2,1,outputlow.f) 
  
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
  tx1 <- cbind(rho*tau2hat[1]*covar.sep(x, X1, d=1/b[1:8], g=g), 
               rho^2*tau2hat[1]*covar.sep(x, X2, d=1/b[1:8], g=g) + tau2hat[2]*covar.sep(x, X2, d=1/b[9:16], g=g))
  
  V1 <- tau2hat[1]*covar.sep(X1, d=1/b[1:8], g=g)
  V12 <- rho*tau2hat[1]*covar.sep(X1, X2, d=1/b[1:8], g=0)
  V2 <- rho^2*tau2hat[1]*covar.sep(X2, d=1/b[1:8], g=g) + tau2hat[2]*covar.sep(X2, d=1/b[9:16], g=g)
  
  V_2 <- rbind(cbind(V1, V12), cbind(t(V12), V2))
  
  mx1 <- tx1 %*% solve(V_2) %*% c(Y1, Y2)
  
  ### posterior variance ###
  koh.var1 <- pmax(0, diag(tau2hat[2]*covar.sep(as.matrix(x), d=1/b[9:16], g=g) + tau2hat[1]*rho^2*covar.sep(as.matrix(x), d=1/b[1:8], g=g) - tx1 %*% solve(V_2)%*%t(tx1)))
  
  return(list(mu=mx1, sig2=koh.var1))
}

costmatk <- list(NA)
rmsematk <- list(NA)
### synthetic function ###
borehole <- function(xx)
{
  rw <- xx[1]
  r  <- xx[2]
  Tu <- xx[3]
  Hu <- xx[4]
  Tl <- xx[5]
  Hl <- xx[6]
  L  <- xx[7]
  Kw <- xx[8]
  
  frac1 <- 2 * pi * Tu * (Hu-Hl)
  
  frac2a <- 2*L*Tu / (log(r/rw)*rw^2*Kw)
  frac2b <- Tu / Tl
  frac2 <- log(r/rw) * (1+frac2a+frac2b)
  
  y <- frac1 / frac2
  return(y)
}

boreholelow <- function(xx)
{ 
  rw <- xx[1]
  r  <- xx[2]
  Tu <- xx[3]
  Hu <- xx[4]
  Tl <- xx[5]
  Hl <- xx[6]
  L  <- xx[7]
  Kw <- xx[8]
  
  frac1 <- 5 * Tu * (Hu-Hl) #+ (Tu * Kw) # Tu * Kw is added
  
  frac2a <- 2*L*Tu / (log(r/rw)*rw^2*Kw)
  frac2b <- Tu / Tl
  frac2 <- log(r/rw) * (1.5+frac2a+frac2b)
  
  y <- frac1 / frac2
  return(y)
}

output.f <- function(x){
  factor_range <- list("rw" = c(0.05, 0.15), "r" = c(100, 50000),
                       "Tu" = c(63070, 115600), "Hu" = c(990, 1110),
                       "Tl" = c(63.1, 116), "Hl" = c(700, 820),
                       "L" = c(1120, 1680), "Kw" = c(9855, 12045))
  for(i in 1:length(factor_range)) x[i] <- factor_range[[i]][1] + x[i] * diff(factor_range[[i]])
  borehole(x[1:8])
} 

outputlow.f <- function(x){
  factor_range <- list("rw" = c(0.05, 0.15), "r" = c(100, 50000),
                       "Tu" = c(63070, 115600), "Hu" = c(990, 1110),
                       "Tl" = c(63.1, 116), "Hl" = c(700, 820),
                       "L" = c(1120, 1680), "Kw" = c(9855, 12045))
  for(i in 1:length(factor_range)) x[i] <- factor_range[[i]][1] + x[i] * diff(factor_range[[i]])
  boreholelow(x[1:8])
} 

### training data ###
n1 <- 40; n2 <- 20
d <- 8
for(kk in 1:10){
  set.seed(kk)
  X1 <- maximinLHS(n1, d)
  X2 <- maximinLHS(n2, d)
  
  NestDesign <- NestedDesignBuild(design = list(X1,X2))
  
  X1 <- NestDesign$PX
  X2 <- ExtractNestDesign(NestDesign,2)
  
  y1 <- apply(X1,1,outputlow.f)
  y2 <- apply(X2,1,output.f)
  
  ### test data ###
  x <- maximinLHS(100, d)
  
  ### KOH method ###
  fit.KOH2 <- fit.KOH(X1, X2, y1, y2)
  pred.KOH2 <- pred.KOH(fit.KOH2, x)
  mx1 <- pred.KOH2$mu
  koh.var1 <- pred.KOH2$sig2
  
  sqrt(mean((mx1-apply(x,1,output.f))^2)) # KOH
  
  
  ### IMSPE ###
  Icurrent <- mean(koh.var1) # current IMSPE
  Icurrent
  
  ### Add 1 points and calculate IMSPE ###
  IcandKOH1 <- c(rep(0, nrow(x))) # IMSPE candidates
  IcandKOH2 <- c(rep(0, nrow(x))) # IMSPE candidates
  
  for(i in 1:length(IcandKOH1)){ # no true, no need to fit just pred
    IcandKOH1[i] <- IMSPEKOH1(x, x[i,], fit.KOH2, level=1)
  }
  for(i in 1:length(IcandKOH2)){ # no true, no need to fit just pred
    IcandKOH2[i] <- IMSPEKOH1(x, x[i,], fit.KOH2, level=2)
  }
  
  which.min(IcandKOH1)
  which.min(IcandKOH2)
  
  ### Fast update; Equation 6.6. in Surrogates ###
  ### ALC; How much can be improved. Equation 6.6. in Surrogates ###
  alcfast <- c(Icurrent - IcandKOH1[which.min(IcandKOH1)], Icurrent - IcandKOH2[which.min(IcandKOH2)])
  alcfast
  
  ### cost; 1, 2, 3 ###
  which.max(alcfast/c(5,10))
  alcfast/c(5,10)
  
  
  chosen <- matrix(0, ncol=2)
  chosen[1,1] <- which.max(alcfast/c(5,10))
  chosen[1,2] <- which.min(cbind(IcandKOH1, IcandKOH2)[,chosen[1,1]])
  
  
  nonlinear.cost <- 0
  nonlinear.error <- sqrt(mean((mx1-apply(x,1,output.f))^2))
  
  Iselect <- IMSPEKOHselect1(x, matrix(x[chosen[nrow(chosen),2],], nrow=1), fit.KOH2, level=chosen[nrow(chosen),1]) 
  
  
  #################
  ### Add point ###
  #################
  while(nonlinear.cost[length(nonlinear.cost)] < 100){ # if total cost is less than the budget
    
    ### predictive ###
    pred.KOH2 <- pred.KOH(Iselect, x)
    mx1 <- pred.KOH2$mu
    koh.var1 <- pred.KOH2$sig2
    
    ### RMSE ###  
    nonlinear.error <- c(nonlinear.error, sqrt(mean((mx1-apply(x,1,output.f))^2))) # closed form
    if(chosen[nrow(chosen),1] == 1){
      nonlinear.cost[length(nonlinear.cost)+1] <- nonlinear.cost[length(nonlinear.cost)]+5
    }else{
      nonlinear.cost[length(nonlinear.cost)+1] <- nonlinear.cost[length(nonlinear.cost)]+10
    }
    
    #############
    ### IMSPE ###
    #############
    Icurrent <- mean(koh.var1)
    
    ### Add 1 points and calculate IMSPE ###
    IcandKOH1 <- c(rep(0, nrow(x))) # IMSPE candidates
    IcandKOH2 <- c(rep(0, nrow(x))) # IMSPE candidates
    
    for(i in 1:length(IcandKOH1)){ # no true, no need to fit just pred
      if(any(chosen[,2]==i)){IcandKOH1[i] <- 0}else{
        IcandKOH1[i] <- IMSPEKOH1(x, x[i,], Iselect, level=1)
      }
    }
    for(i in 1:length(IcandKOH2)){ # no true, no need to fit just pred
      if(any(chosen[,2]==i)){IcandKOH2[i] <- 0}else{
        IcandKOH2[i] <- IMSPEKOH1(x, x[i,], Iselect, level=2)
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
    which.max(alcfast/c(5,10))
    alcfast/c(5,10)
    
    
    chosen <- rbind(chosen, c(which.max(alcfast/c(5,10)), which.min(cbind(IcandKOH1, IcandKOH2)[,which.max(alcfast/c(5,10))])))
    Iselect <- IMSPEKOHselect1(x, matrix(x[chosen[nrow(chosen),2],], nrow=1), Iselect, level=chosen[nrow(chosen),1])
    
    if(nonlinear.cost[length(nonlinear.cost)] >=100){break}
    
  }
  
  
  ### Save results ###
  costmatk[[kk]] <- nonlinear.cost
  rmsematk[[kk]] <- nonlinear.error
}
costmatk
rmsematk



