source('/Users/junoh/Downloads/code_technometrics/Internal_functions.R')
source('/Users/junoh/Downloads/code_technometrics/Sequential_algorithms.R')
library('DiceKriging')
library('MuFiCokriging')
library('lhs')
library('rgenoud')
library('mcmc')

#-----------------------------------------------------------------------------#
#----------------#             sequential co-kriging         #----------------#
#-----------------------------------------------------------------------------#

costmatco <- list(NA)
rmsematco <- list(NA)

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

eps <- sqrt(.Machine$double.eps)

#--- Initial experimental design sets
# iter <- c(12,18,10,18,10,11,13,12,25,10)
for(kk in 1:10){
  set.seed(kk) 
  
  n1 <- 15; n2 <- 8
  d <- 4
  
  X1 <- maximinLHS(n1, d)
  X2 <- maximinLHS(n2, d)
  
  NestDesign <- NestedDesignBuild(design = list(X1,X2))
  
  X1 <- NestDesign$PX
  X2 <- ExtractNestDesign(NestDesign,2)
  
  y1 <- apply(X1,1,park91alc)
  y2 <- apply(X2,1,park91a)
  
  x <- maximinLHS(100, d)
  
  #--- Multi-fidelity co-kriging model building
  mymodel <- MuFicokm(
    formula = list(~1,~1), 
    MuFidesign = NestDesign, 
    response = list(y1,y2), 
    lower=eps, upper=0.1,
    coef.trend = list(0,c(0,0)),
    nlevel = 2,
    covtype = "gauss"
  )
  res <- summary(mymodel)

  predictions <- predict(
    object = mymodel, 
    newdata = x,
    type="SK")
  ###
  rmsematco[[kk]] <- sqrt(mean((predictions$mean-apply(x,1,park91a))^2)) # Cokm
  
  ## One point at-a-time sequential cokriging (see Section 3.1). ##
  B <- 10		#-B: ratio of computational costs between level 1 and 2
  
  for(i in 1:35){ #iter[kk]
    ###
    set.seed(kk) ### 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
    
    n1 <- 15; n2 <- 8
    
    X1 <- maximinLHS(n1, d)
    X2 <- maximinLHS(n2, d)
    
    NestDesign <- NestedDesignBuild(design = list(X1,X2))
    
    X1 <- NestDesign$PX
    X2 <- ExtractNestDesign(NestDesign,2)
    
    y1 <- apply(X1,1,park91alc)
    y2 <- apply(X2,1,park91a)
    
    #--- Multi-fidelity co-kriging model building
    mymodel <- MuFicokm(
      formula = list(~1,~1), 
      MuFidesign = NestDesign, 
      response = list(y1,y2), 
      coef.trend = list(0,c(0,0)),
      lower=eps, upper=0.1,
      nlevel = 2,
      covtype = "gauss"
    )
    res <- summary(mymodel)

    predictions <- predict(object = mymodel, newdata = x, type="SK")
    
    ###
    niter <- i ### 12,18,10,18,10,11,13,12,25,10
    
    cokm_varmax <- one_step_cokm_varmax(
      model = mymodel,
      B = B,
      xpred = x,
      yreal = park91a(x),
      myfunctions = list(park91alc, park91a),
      niter = niter,
      param.estim = TRUE,
      error.compute = FALSE,
      error.LOO = FALSE,
      ponderation = FALSE)
    
    ###
    rmsematco[[kk]][i+1] <- sqrt(mean((cokm_varmax$ypredseq$mean-apply(x,1,park91a))^2)) # Cokm
  }
  
  # cokm_varmax$CoutSave <- 2*cokm_varmax$CoutSave
  ###
  costmatco[[kk]] <- c(0,cokm_varmax$CoutSave)
}
costmatco
rmsematco

