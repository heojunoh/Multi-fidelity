library(lhs)
library(laGP)

costmatc <- list(NA)
rmsematc <- list(NA)
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
mcsam <- c(4,4,4,4,4,4,4,4,4,4)
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
  
  
  ### closed ###
  fit.closed <- RNAmf(X1, y1, X2, y2, kernel="sqex", constant=TRUE)
  predy <- predRNAmf(fit.closed, x)$mu
  predsig2 <- predRNAmf(fit.closed, x)$sig2
  
  ### RMSE ###
  sqrt(mean((predy-apply(x,1,output.f))^2)) # closed form
  
  
  ### IMSPE ###
  Icurrent <- mean(predsig2) # current IMSPE
  Icurrent
  
  ### Add 1 points and calculate IMSPE ###
  intgvr <- integvar(x, fit.closed, mc.sample=mcsam[kk])
  
  Icand1fast <- intgvr$intvar1
  Icand2fast <- intgvr$intvar2
  
  which.min(Icand1fast)
  which.min(Icand2fast)
  
  ### Fast update; Equation 6.6. in Surrogates ###
  ### ALC; How much can be improved. Equation 6.6. in Surrogates ###
  alcfast <- c(Icurrent - Icand1fast[which.min(Icand1fast)], Icurrent - Icand2fast[which.min(Icand2fast)])
  alcfast
  
  ### cost; 1, 2, 3 ###
  which.max(alcfast/c(5,(5+10)))
  alcfast/c(5,(5+10))
  
  
  chosen <- matrix(0, ncol=2)
  chosen[1,1] <- which.max(alcfast/c(5,(5+10)))
  chosen[1,2] <- which.min(cbind(Icand1fast, Icand2fast)[,chosen[1,1]])
  
  
  nonlinear.cost <- 0
  nonlinear.error <- sqrt(mean((predy-apply(x,1,output.f))^2))
  
  Iselect <- IMSPEselect1(x[chosen[nrow(chosen),2],], fit.closed, level=chosen[nrow(chosen),1])
  
  
  
  #################
  ### Add point ###
  #################
  while(nonlinear.cost[length(nonlinear.cost)] < 100){ # if total cost is less than the budget
    
    ### closed ###
    predy <- predRNAmf(Iselect$fit, x)$mu
    predsig2 <- predRNAmf(Iselect$fit, x)$sig2
    
    ### RMSE ###
    nonlinear.error <- c(nonlinear.error, sqrt(mean((predy-apply(x,1,output.f))^2))) # closed form
    if(chosen[nrow(chosen),1] == 1){
      nonlinear.cost[length(nonlinear.cost)+1] <- nonlinear.cost[length(nonlinear.cost)]+5
    }else{
      nonlinear.cost[length(nonlinear.cost)+1] <- nonlinear.cost[length(nonlinear.cost)]+(5+10)
    }
    
    #############
    ### IMSPE ###
    #############
    Icurrent <- mean(predsig2)
    
    ### Add 1 points to the low-fidelity data ###
    intgvr <- integvar(x, Iselect$fit, mc.sample=mcsam[kk])
    
    Icand1fast <- intgvr$intvar1
    Icand2fast <- intgvr$intvar2
    
    for(i in 1:nrow(x)){ # no true, no need to fit just pred
      if(any(chosen[,2]==i)){
        Icand1fast[i] <- max(Icand1fast)
        Icand2fast[i] <- max(Icand2fast)
      }
    }
    
    which.min(Icand1fast)
    which.min(Icand2fast)
    
    ### Fast update; Equation 6.6. in Surrogates ###
    ### ALC; How much can be improved. Equation 6.6. in Surrogates ###
    alcfast <- c(Icurrent - Icand1fast[which.min(Icand1fast)], Icurrent - Icand2fast[which.min(Icand2fast)])
    alcfast
    
    ### cost; 1, 2, 3 ###
    which.max(alcfast/c(5,(5+10)))
    alcfast/c(5,(5+10))
    
    
    chosen <- rbind(chosen, c(which.max(alcfast/c(5,(5+10))), which.min(cbind(Icand1fast, Icand2fast)[,which.max(alcfast/c(5,(5+10)))])))
    Iselect <- IMSPEselect1(x[chosen[nrow(chosen),2],], Iselect$fit, level=chosen[nrow(chosen),1])
    
    if(nonlinear.cost[length(nonlinear.cost)] >= 100){break}
    
  }
  
  
  ### Save results ###
  costmatc[[kk]] <- nonlinear.cost
  rmsematc[[kk]] <- nonlinear.error
}
costmatc
rmsematc



