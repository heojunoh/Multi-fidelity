library(lhs)
library(laGP)

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

  ### test data ###
  x <- maximinLHS(100, d)


  ### closed ###
  fit.closed <- RNAmf(X1, y1, X2, y2, kernel="sqex", constant=TRUE)
  predy <- predRNAmf(fit.closed, x)$mu
  predsig2 <- predRNAmf(fit.closed, x)$sig2

  ### RMSE ###
  sqrt(mean((predy-apply(x,1,park91a))^2)) # closed form


  ### IMSPE ###
  Icurrent <- mean(predsig2) # current IMSPE
  Icurrent

  ### Add 1 points and calculate IMSPE ###
  which.max(predsig2)

  intgvr <- integvar(matrix(x[which.max(predsig2),], nrow=1), fit.closed, mc.sample=10)

  Icand1fast <- intgvr$intvar1
  Icand2fast <- intgvr$intvar2

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
  nonlinear.error <- sqrt(mean((predy-apply(x,1,park91a))^2))

  Iselect <- IMSPEselect1(x[chosen[nrow(chosen),2],], fit.closed, level=chosen[nrow(chosen),1])



  #################
  ### Add point ###
  #################
  while(nonlinear.cost[length(nonlinear.cost)] < 100){ # if total cost is less than the budget

    ### closed ###
    predy <- predRNAmf(Iselect$fit, x)$mu
    predsig2 <- predRNAmf(Iselect$fit, x)$sig2

    ### RMSE ###
    nonlinear.error <- c(nonlinear.error, sqrt(mean((predy-apply(x,1,park91a))^2))) # closed form
    if(chosen[nrow(chosen),1] == 1){
      nonlinear.cost[length(nonlinear.cost)+1] <- nonlinear.cost[length(nonlinear.cost)]+1
    }else{
      nonlinear.cost[length(nonlinear.cost)+1] <- nonlinear.cost[length(nonlinear.cost)]+(1+10)
    }

    #############
    ### IMSPE ###
    #############
    Icurrent <- mean(predsig2)

    ### Add 1 points to the low-fidelity data ###
    intgvr <- integvar(matrix(x[which.max(predsig2),], nrow=1), Iselect$fit, mc.sample=10)

    Icand1fast <- intgvr$intvar1
    Icand2fast <- intgvr$intvar2

    ### Fast update; Equation 6.6. in Surrogates ###
    ### ALC; How much can be improved. Equation 6.6. in Surrogates ###
    alcfast <- c(Icurrent - Icand1fast, Icurrent - Icand2fast)
    alcfast

    ### cost; 1, 2, 3 ###
    which.max(alcfast/c(1,(1+10)))
    alcfast/c(1,(1+10))


    chosen <- rbind(chosen, c(which.max(alcfast/c(1,(1+10))), which.max(predsig2)))
    Iselect <- IMSPEselect1(x[chosen[nrow(chosen),2],], Iselect$fit, level=chosen[nrow(chosen),1])

    if(nonlinear.cost[length(nonlinear.cost)] >= 100){break}

  }


  ### Save results ###
  costmatc2[[kk]] <- nonlinear.cost
  rmsematc2[[kk]] <- nonlinear.error
}
costmatc2
rmsematc2



