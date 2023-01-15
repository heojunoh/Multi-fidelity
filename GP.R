library(plgp)
eps <- sqrt(.Machine$double.eps) 

GP <- function(X, y, g=eps, center=TRUE){
  
  n <- length(y)
  # if(center) y <- scale(y, scale=FALSE)
  
  nlsep <- function(par, X, Y) 
  {
    n <- length(Y)
    theta <- par # lengthscale
    K <- covar.sep(X, d=theta, g=g)
    Ki <- solve(K+diag(g,n))
    ldetK <- determinant(K, logarithm=TRUE)$modulus
    ll <- - (n/2)*log(t(Y) %*% Ki %*% Y) - (1/2)*ldetK
    return(drop(-ll))
  }
  
  outg <- optim(c(rep(0.1, ncol(X))), nlsep, method="L-BFGS-B", 
                lower=c(rep(0.001*sqrt(ncol(X)), ncol(X))), upper=c(rep(1000*sqrt(ncol(X)), ncol(X))), X=X, Y=y) 
  
  K <- covar.sep(X, d=outg$par, g=g)
  Ki <- solve(K+diag(g,n))
  tau2hat <- drop(t(y) %*% Ki %*% y / n)
  
  return(list(theta = outg$par, g=g, Ki=Ki, X = X, y = y, tau2hat=tau2hat))
}

pred.GP <- function(fit, xnew){
  
  xnew <- as.matrix(xnew)
  
  Ki <- fit$Ki
  theta <- fit$theta
  nu <- fit$nu
  theta <- fit$theta
  g <- fit$g
  X <- fit$X
  y <- fit$y
  tau2hat <- fit$tau2hat
  
  KXX <- covar.sep(as.matrix(xnew), d=theta, g=g) 
  KX <- covar.sep(xnew, X, d=theta, g=0)
  mup2 <- KX %*% Ki %*% y
  Sigmap2 <- pmax(0, diag(tau2hat*(KXX - KX %*% Ki %*% t(KX))))
 
  return(list(mu=mup2, sig2=Sigmap2))
}
