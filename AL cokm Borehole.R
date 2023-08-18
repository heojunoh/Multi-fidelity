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

eps <- sqrt(.Machine$double.eps)

#--- Initial experimental design sets
# iter <- c(12,18,10,18,10,11,13,12,25,10)
for(kk in 1:10){
  set.seed(kk) 
  
  n1 <- 40; n2 <- 20
  d <- 8
  
  X1 <- maximinLHS(n1, d)
  X2 <- maximinLHS(n2, d)
  
  NestDesign <- NestedDesignBuild(design = list(X1,X2))
  
  X1 <- NestDesign$PX
  X2 <- ExtractNestDesign(NestDesign,2)
  
  y1 <- apply(X1,1,outputlow.f)
  y2 <- apply(X2,1,output.f)
  
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
  rmsematco[[kk]] <- sqrt(mean((predictions$mean-apply(x,1,output.f))^2)) # Cokm
  
  ## One point at-a-time sequential cokriging (see Section 3.1). ##
  B <- 2		#-B: ratio of computational costs between level 1 and 2
  
  for(i in 1:20){ #iter[kk]
    ###
    set.seed(kk) ### 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
    
    n1 <- 40; n2 <- 20
    
    X1 <- maximinLHS(n1, d)
    X2 <- maximinLHS(n2, d)
    
    NestDesign <- NestedDesignBuild(design = list(X1,X2))
    
    X1 <- NestDesign$PX
    X2 <- ExtractNestDesign(NestDesign,2)
    
    y1 <- apply(X1,1,outputlow.f)
    y2 <- apply(X2,1,output.f)
    
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
      yreal = apply(x,1,output.f),
      # myfunctions = list(park91alc, park91a),
      niter = niter,
      param.estim = TRUE,
      error.compute = FALSE,
      error.LOO = FALSE,
      ponderation = FALSE)
    
    ###
    rmsematco[[kk]][i+1] <- sqrt(mean((cokm_varmax$ypredseq$mean-apply(x,1,output.f))^2)) # Cokm
  }
  
  # cokm_varmax$CoutSave <- 2*cokm_varmax$CoutSave
  ###
  costmatco[[kk]] <- c(0,5*cokm_varmax$CoutSave)
}

costmatco
rmsematco

