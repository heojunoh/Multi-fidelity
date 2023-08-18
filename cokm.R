##########################
### Le Gartiet package ###
##########################

# install.packages("devtools")
library(devtools)
# install_github("cran/MuFiCokriging")
library(MuFiCokriging)
library(lhs)
library(plgp)


IMSPEcokm <- function(x, newx, fit){
  sumcokm <- summary(fit)
  
  k2 <- k3 <- rep(0, length(x))
  
  # k_n_1, no h and H because no l-1
  f1 <- 1 
  F1 <- matrix(f1, ncol=1, nrow=nrow(fit$cok[[1]]@X))
  rr1 <- covar.sep(matrix(x), d=(sumcokm$Cov.Val[[1]])^2, g=eps)
  r1 <- covar.sep(x, rbind(fit$cok[[1]]@X, newx), d=(sumcokm$Cov.Val[[1]])^2, g=eps)
  R1 <- covar.sep(rbind(fit$cok[[1]]@X, newx), d=(sumcokm$Cov.Val[[1]])^2, g=eps)
  
  k1 <- diag(sumcokm$Var.Val[[1]]*( rr1 - r1 %*% solve(R1) %*% t(r1) ))
  
  
  # k_n_2
  # F and f are 1; constant trend, beta is in sum.muficokm$Trend.Val
  # H, R, F
  f2 <- 1 
  F2 <- matrix(f2, ncol=1, nrow=nrow(fit$cok[[2]]@X))
  H2 <- rbind(fit$cok[[2]]@AR.F, predict(fit, newx, "SK")$mux[[1]]) %*% f2 # y1(D2 & newx)
  h2 <- matrix(predict(fit, x, "SK")$mux[[1]] * f2)   # predictive mean at level 1 * f2
  rr2 <- covar.sep(matrix(x), d=(sumcokm$Cov.Val[[2]])^2, g=eps)
  r2 <- covar.sep(x, rbind(fit$cok[[2]]@X, newx), d=(sumcokm$Cov.Val[[2]])^2, g=eps)
  R2 <- covar.sep(rbind(fit$cok[[2]]@X, newx), d=(sumcokm$Cov.Val[[2]])^2, g=eps)
  
  sig2rho1 <- sumcokm$Rho.Val[[1]]^2 + solve(t(H2) %*% solve(R2) %*% H2)
  sig2delta2 <- sumcokm$Var.Val[[2]]*diag( rr2 - cbind(h2, r2) %*% solve( rbind(cbind(0, t(H2)), cbind(H2, R2)) ) %*% t(cbind(h2, r2)) )
  
  k2 <- k1 * drop(sig2rho1) + sig2delta2
  
  
  # k_n_3
  # F and f are 1; constant trend, beta is in sum.muficokm$Trend.Val
  # H, R, F
  f3 <- 1 
  F3 <- matrix(f3, ncol=1, nrow=nrow(fit$cok[[3]]@X))
  H3 <- rbind(fit$cok[[3]]@AR.F, predict(fit, newx, "SK")$mux[[2]]) %*% f3 # y2(D3 & newx)
  h3 <- matrix(predict(fit, x, "SK")$mux[[2]] * f3)   # predictive mean at level 2 * f3
  rr3 <- covar.sep(matrix(x), d=(sumcokm$Cov.Val[[3]])^2, g=eps)
  r3 <- covar.sep(x, rbind(fit$cok[[3]]@X, newx), d=(sumcokm$Cov.Val[[3]])^2, g=eps)
  R3 <- covar.sep(rbind(fit$cok[[3]]@X, newx), d=(sumcokm$Cov.Val[[3]])^2, g=eps)
  
  sig2rho2 <- sumcokm$Rho.Val[[2]]^2 + solve(t(H3) %*% solve(R3) %*% H3)
  sig2delta3 <- sumcokm$Var.Val[[3]]*diag( rr3 - cbind(h3, r3) %*% solve( rbind(cbind(0, t(H3)), cbind(H3, R3)) ) %*% t(cbind(h3, r3)) )
  
  k3 <- k2 * drop(sig2rho2) + sig2delta3
  
  return(list(k1=mean(k1), k2=mean(k2), k3=mean(k3)))
}


### This is updating model everytime.
IMSPEcokmselect <- function(x, newx, fit, level){
  
  ### choose level 1 ###
  if(level == 1){
    y1.select <- fl(newx, l=1)
    fit$Dnest$PX <- rbind(fit$Dnest$PX, newx)
    
    fit.new <- MuFicokm(formula = list(~1,~1,~1), MuFidesign = fit$Dnest, covtype="gauss",
                        # coef.trend = list(0,c(0,0),c(0,0)),
                        lower=eps, upper=0.1,
                        response = list(c(fit$cok[[1]]@y, y1.select),fit$cok[[2]]@y,fit$cok[[3]]@y), nlevel = 3)
  }
  
  ### choose level 2 ###
  if(level == 2){
    y1.select <- fl(newx, l=1)
    y2.select <- fl(newx, l=3)
    fit$Dnest$PX <- rbind(fit$Dnest$PX, newx)
    fit$Dnest$ind[[1]] <- c(fit$Dnest$ind[[1]], nrow(fit$Dnest$PX))
    fit$Dnest$n[1] <- fit$Dnest$n[1]+1
    
    fit.new <- MuFicokm(formula = list(~1,~1,~1), MuFidesign = fit$Dnest, covtype="gauss",
                        # coef.trend = list(0,c(0,0),c(0,0)),
                        lower=eps, upper=0.1,
                        response = list(c(fit$cok[[1]]@y, y1.select),c(fit$cok[[2]]@y, y2.select),fit$cok[[3]]@y), nlevel = 3)
  }
  
  ### choose level 3 ###
  if(level == 3){
    y1.select <- fl(newx, l=1)
    y2.select <- fl(newx, l=3)
    y3.select <- fl(newx, l=5)
    fit$Dnest$PX <- rbind(fit$Dnest$PX, newx)
    fit$Dnest$ind[[1]] <- c(fit$Dnest$ind[[1]], nrow(fit$Dnest$PX))
    fit$Dnest$ind[[2]] <- c(fit$Dnest$ind[[2]], length(fit$Dnest$ind[[1]]))
    fit$Dnest$n[1] <- fit$Dnest$n[1]+1
    fit$Dnest$n[2] <- fit$Dnest$n[2]+1
    
    fit.new <- MuFicokm(formula = list(~1,~1,~1), MuFidesign = fit$Dnest, covtype="gauss",
                        # coef.trend = list(0,c(0,0),c(0,0)),
                        lower=eps, upper=0.1,
                        response = list(c(fit$cok[[1]]@y, y1.select),c(fit$cok[[2]]@y, y2.select),c(fit$cok[[3]]@y, y3.select)), nlevel = 3)
  }
  return(list(fit=fit.new))
}


IMSPEcokm1 <- function(x, newx, fit){
  sumcokm <- summary(fit)
  
  k2 <- rep(0, length(x))
  
  # k_n_1, no h and H because no l-1
  f1 <- 1 
  F1 <- matrix(f1, ncol=1, nrow=nrow(fit$cok[[1]]@X))
  rr1 <- covar.sep(matrix(x), d=(sumcokm$Cov.Val[[1]])^2, g=eps)
  r1 <- covar.sep(x, rbind(fit$cok[[1]]@X, newx), d=(sumcokm$Cov.Val[[1]])^2, g=eps)
  R1 <- covar.sep(rbind(fit$cok[[1]]@X, newx), d=(sumcokm$Cov.Val[[1]])^2, g=eps)
  
  k1 <- diag(sumcokm$Var.Val[[1]]*( rr1 - r1 %*% solve(R1) %*% t(r1) ))
  
  # k_n_2
  # F and f are 1; constant trend, beta is in sum.muficokm$Trend.Val
  # H, R, F
  f2 <- 1 
  F2 <- matrix(f2, ncol=1, nrow=nrow(fit$cok[[2]]@X))
  H2 <- rbind(fit$cok[[2]]@AR.F, predict(fit, newx, "SK")$mux[[1]]) %*% f2 # y1(D2 & newx)
  h2 <- matrix(predict(fit, x, "SK")$mux[[1]] * f2)   # predictive mean at level 1 * f2
  rr2 <- covar.sep(matrix(x), d=(sumcokm$Cov.Val[[2]])^2, g=eps)
  r2 <- covar.sep(x, rbind(fit$cok[[2]]@X, newx), d=(sumcokm$Cov.Val[[2]])^2, g=eps)
  R2 <- covar.sep(rbind(fit$cok[[2]]@X, newx), d=(sumcokm$Cov.Val[[2]])^2, g=eps)
  
  sig2rho1 <- sumcokm$Rho.Val[[1]]^2 + solve(t(H2) %*% solve(R2) %*% H2)
  sig2delta2 <- sumcokm$Var.Val[[2]]*diag( rr2 - cbind(h2, r2) %*% solve( rbind(cbind(0, t(H2)), cbind(H2, R2)) ) %*% t(cbind(h2, r2)) )
  
  k2 <- k1 * drop(sig2rho1) + sig2delta2
  
  return(list(k1=mean(k1), k2=mean(k2)))
}


### This is updating model everytime.
IMSPEcokmselect1 <- function(x, newx, fit, level){
  
  ### choose level 1 ###
  if(level == 1){
    y1.select <- f1(newx)
    fit$Dnest$PX <- rbind(fit$Dnest$PX, newx)
    
    fit.new <- MuFicokm(formula = list(~1,~1), MuFidesign = fit$Dnest, covtype="gauss",
                        # coef.trend = list(0,c(0,0),c(0,0)),
                        lower=eps, upper=1/nrow(fit$Dnest$PX),
                        response = list(c(fit$cok[[1]]@y, y1.select),fit$cok[[2]]@y), nlevel = 2)
  }
  
  ### choose level 2 ###
  if(level == 2){
    y1.select <- f1(newx)
    y2.select <- f2(newx)
    fit$Dnest$PX <- rbind(fit$Dnest$PX, newx)
    fit$Dnest$ind[[1]] <- c(fit$Dnest$ind[[1]], nrow(fit$Dnest$PX))
    fit$Dnest$n[1] <- fit$Dnest$n[1]+1
    
    fit.new <- MuFicokm(formula = list(~1,~1), MuFidesign = fit$Dnest, covtype="gauss",
                        # coef.trend = list(0,c(0,0),c(0,0)),
                        lower=eps, upper=1/nrow(fit$Dnest$PX),
                        response = list(c(fit$cok[[1]]@y, y1.select),c(fit$cok[[2]]@y, y2.select)), nlevel = 2)
  }
  
  return(fit=fit.new)
}

