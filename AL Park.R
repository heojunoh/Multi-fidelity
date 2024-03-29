library(lhs)
library(laGP)
library(plgp)
library(MuFiCokriging)
library(doParallel)
library(foreach)
library(RNAmf)

eps <- sqrt(.Machine$double.eps)
crps <- function(x, mu, sig2){ # The smaller, the better (0 to infinity)
  if(any(sig2==0)) sig2[sig2==0] <- eps
  -sqrt(sig2)*(1/sqrt(pi)-2*dnorm((x-mu)/sqrt(sig2))-(x-mu)/sqrt(sig2)*(2*pnorm((x-mu)/sqrt(sig2))-1))
}

costmatc <- list(NA)
rmsematc <- list(NA)
crpsmatc <- list(NA)
time.each <- rep(0,10)
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
n1 <- 40; n2 <- 20
d <- 4

for(kk in 1:10){
  time.start <- proc.time()[3]
  set.seed(kk)
  print(kk)
  X1 <- maximinLHS(n1, d)
  X2 <- maximinLHS(n2, d)

  NestDesign <- NestedDesignBuild(design = list(X1,X2))

  X1 <- NestDesign$PX
  X2 <- ExtractNestDesign(NestDesign,2)

  y1 <- apply(X1,1,park91alc)
  y2 <- apply(X2,1,park91a)

  ### test data ###
  x <- maximinLHS(1000, d)


  ### closed ###
  fit.closed <- RNAmf(X1, y1, X2, y2, kernel="sqex", constant=TRUE)
  predy <- predRNAmf(fit.closed, x)$mu
  predsig2 <- predRNAmf(fit.closed, x)$sig2

  ### RMSE ###
  sqrt(mean((predy-apply(x,1,park91a))^2)) # closed form
  
  park.cost <- 0
  park.error <- sqrt(mean((predy-apply(x,1,park91a))^2))
  park.crps <- mean(crps(apply(x,1,park91a), predy, predsig2))

  Iselect <- ALM_two_level(fit.closed, c(1,6), list(park91alc, park91a), parallel=TRUE, ncore=10)


  #################
  ### Add point ###
  #################
  while(park.cost[length(park.cost)] < 50){ # if total cost is less than the budget
    ### closed ###
    predy <- predRNAmf(Iselect$fit, x)$mu
    predsig2 <- predRNAmf(Iselect$fit, x)$sig2

    ### RMSE ###
    park.error <- c(park.error, sqrt(mean((predy-apply(x,1,park91a))^2))) # RMSE
    park.crps <- c(park.crps, mean(crps(apply(x,1,park91a), predy, predsig2))) # CRPS
    if(Iselect$chosen$level == 1){
      park.cost[length(park.cost)+1] <- park.cost[length(park.cost)]+1
    }else{
      park.cost[length(park.cost)+1] <- park.cost[length(park.cost)]+(1+6)
    }
    print(park.cost[length(park.cost)])
    print(park.error[length(park.error)])
    
    if(park.cost[length(park.cost)] >= 50){break}
    
    ### update the next point ###
    Iselect <- ALM_two_level(Iselect$fit, c(1,6), list(park91alc, park91a), parallel=TRUE, ncore=10)
    # save.image("C:/Users/heojunoh/Desktop/Park AL 1,6.RData")
  }


  ### Save results ###
  costmatc[[kk]] <- park.cost
  rmsematc[[kk]] <- park.error
  crpsmatc[[kk]] <- park.crps
  # save.image("C:/Users/heojunoh/Desktop/Park AL 1,6.RData")
  
  time.each[kk] <- proc.time()[3]- time.start
}
costmatc
rmsematc
crpsmatc
time.each



