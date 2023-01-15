library(plgp)
source("GP.R")
eps <- sqrt(.Machine$double.eps) 

KOH <- function(X, Y2, Y1, g=eps, center=TRUE){
  
  n <- length(Y2)
  # if(center) y <- scale(y, scale=FALSE)
  
  nlsep2 <- function(par, X, Y2, Y1) 
  {
    n <- length(Y2)
    theta <- par # lengthscale
    Y <- Y2-theta[ncol(X)+1]*Y1
    K <- covar.sep(X, d=theta[1:ncol(X)], g=g)
    Ki <- solve(K+diag(g,n))
    ldetK <- determinant(K, logarithm=TRUE)$modulus
    ll <- - (n/2)*log(t(Y) %*% Ki %*% Y) - (1/2)*ldetK
    return(drop(-ll))
  }
  
  outg <- optim(c(rep(0.1, ncol(X)), 0.1), # For lengthscale and rho #
                nlsep2, method="L-BFGS-B", 
                lower=c(rep(0.001*sqrt(ncol(X)), ncol(X)+1)), 
                upper=c(rep(1000*sqrt(ncol(X)), ncol(X)+1)), 
                X=X, Y2=Y2, Y1=Y1) 
  
  K <- covar.sep(X, d=outg$par[1:ncol(X)], g=g)
  Ki <- solve(K+diag(g,n))
  rho <- outg$par[ncol(X)+1]
  y <- Y2-rho*Y1
  tau2hat <- drop(t(y) %*% Ki %*% y / n)
  
  return(list(theta = outg$par[1:ncol(X)], rho = rho, g=g, Ki=Ki, X = X, y = y, tau2hat=tau2hat))
}

