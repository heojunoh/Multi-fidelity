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

#--- functions

f1 <- function(x)
{
  sin(8*pi*x)
}

f2 <- function(x)
{ 
  (x-sqrt(2))*(sin(8*pi*x))^2
}
eps <- sqrt(.Machine$double.eps)

#--- Initial experimental design sets
# iter <- c(12,18,10,18,10,11,13,12,25,10)
for(kk in 1:10){
  set.seed(kk) 
  
  n1 <- 12; n2 <- 9
  
  X1 <- maximinLHS(n1, 1)
  X2 <- maximinLHS(n2, 1)
  
  NestDesign <- NestedDesignBuild(design = list(X1,X2))
  
  X1 <- NestDesign$PX
  X2 <- ExtractNestDesign(NestDesign,2)
  
  y1 <- f1(X1)
  y2 <- f2(X2)
  
  #--- Multi-fidelity co-kriging model building
  mymodel <- MuFicokm(
    formula = list(~1,~1), 
    MuFidesign = NestDesign, 
    response = list(y1,y2), 
    lower=eps, upper=0.1,
    nlevel = 2,
    covtype = "gauss"
  )
  res <- summary(mymodel)
  #--- Sequential design
  x <- seq(0,1,0.01)
  predictions <- predict(
    object = mymodel, 
    newdata = x,
    type="SK")
  ###
  rmsematco[[kk]] <- sqrt(mean((predictions$mean-f2(x))^2)) # Cokm
  
  ## One point at-a-time sequential cokriging (see Section 3.1). ##
  B <- 4		#-B: ratio of computational costs between level 1 and 2
  
  for(i in 1:25){ #iter[kk]
    ###
    set.seed(kk) ### 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
    
    n1 <- 12; n2 <- 9
    
    X1 <- maximinLHS(n1, 1)
    X2 <- maximinLHS(n2, 1)
    
    NestDesign <- NestedDesignBuild(design = list(X1,X2))
    
    X1 <- NestDesign$PX
    X2 <- ExtractNestDesign(NestDesign,2)
    
    y1 <- f1(X1)
    y2 <- f2(X2)
    
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
    #--- Sequential design
    x <- seq(0,1,0.01)
    predictions <- predict(object = mymodel, newdata = x, type="SK")
    
    ###
    niter <- i ### 12,18,10,18,10,11,13,12,25,10
    
    cokm_varmax <- one_step_cokm_varmax(
      model = mymodel,
      B = B,
      xpred = x,
      yreal = f2(x),
      myfunctions = list(f1,f2),
      niter = niter,
      param.estim = TRUE,
      error.compute = FALSE,
      error.LOO = FALSE,
      ponderation = FALSE)
    
    ###
    rmsematco[[kk]][i+1] <- sqrt(mean((cokm_varmax$ypredseq$mean-f2(x))^2)) # Cokm
  }
  
  cokm_varmax$CoutSave <- 2*cokm_varmax$CoutSave
  ###
  costmatco[[kk]] <- c(0,cokm_varmax$CoutSave)
}
costmatco
rmsematco

