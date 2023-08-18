library(lhs)
library(laGP)

costmatc3 <- list(NA)
rmsematc3 <- list(NA)
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
  
  
  ### RNAmf ###
  fit.RNAmf <- RNAmf(X1, y1, X2, y2, kernel="sqex", constant=TRUE)
  predy <- predRNAmf(fit.RNAmf, x)$mu
  predsig2 <- predRNAmf(fit.RNAmf, x)$sig2
  
  ### RMSE ###
  sqrt(mean((predy-apply(x,1,output.f))^2)) # RNAmf
  
  
  ### current uncertainty at each level ###
  unc <- c(max(pred.GP(fit.RNAmf$fit1, x)$sig2), max(predsig2))
  uncind <- c(which.max(pred.GP(fit.RNAmf$fit1, x)$sig2), which.max(predsig2))
  
  ### cost; 1, 2, 3 ###
  which.max(unc/c(5,(5+10)))
  unc/c(5,(5+10))
  
  
  chosen <- matrix(0, ncol=2)
  chosen[1,1] <- which.max(unc/c(5,(5+10)))
  chosen[1,2] <- uncind[which.max(unc/c(5,(5+10)))]
  
  
  nonlinear.cost <- 0
  nonlinear.error <- sqrt(mean((predy-apply(x,1,output.f))^2))
  
  Iselect <- IMSPEselect1(x[chosen[nrow(chosen),2],], fit.RNAmf, level=chosen[nrow(chosen),1])
  
  
  
  #################
  ### Add point ###
  #################
  while(nonlinear.cost[length(nonlinear.cost)] < 100){ # if total cost is less than the budget
    
    ### RNAmf ###
    predy <- predRNAmf(Iselect$fit, x)$mu
    predsig2 <- predRNAmf(Iselect$fit, x)$sig2
    
    ### RMSE ###
    nonlinear.error <- c(nonlinear.error, sqrt(mean((predy-apply(x,1,output.f))^2))) # RNAmf form
    if(chosen[nrow(chosen),1] == 1){
      nonlinear.cost[length(nonlinear.cost)+1] <- nonlinear.cost[length(nonlinear.cost)]+5
    }else{
      nonlinear.cost[length(nonlinear.cost)+1] <- nonlinear.cost[length(nonlinear.cost)]+(5+10)
    }
    
    #############
    ### IMSPE ###
    #############
    
    ### current uncertainty at each level ###
    unc <- c(max(pred.GP(Iselect$fit$fit1, x)$sig2), max(predsig2))
    uncind <- c(which.max(pred.GP(Iselect$fit$fit1, x)$sig2), which.max(predsig2))
    
    ### cost; 1, 2, 3 ###
    which.max(unc/c(5,(5+10)))
    unc/c(5,(5+10))
    
    
    chosen <- rbind(chosen, c(which.max(unc/c(5,(5+10))), uncind[which.max(unc/c(5,(5+10)))]))
    Iselect <- IMSPEselect1(x[chosen[nrow(chosen),2],], Iselect$fit, level=chosen[nrow(chosen),1])
    
    if(nonlinear.cost[length(nonlinear.cost)] >= 100){break}
    
  }
  
  
  ### Save results ###
  costmatc3[[kk]] <- nonlinear.cost
  rmsematc3[[kk]] <- nonlinear.error
}
costmatc3
rmsematc3



