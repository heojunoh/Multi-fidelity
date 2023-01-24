library(plgp)
eps <- sqrt(.Machine$double.eps) 

GP <- function(X, y, g=eps, Xscale=TRUE, Yscale=TRUE){
  
  if(is.null(dim(X))) X <- matrix(X, ncol = 1)
  
  if(Xscale){
    X <- scale(X, center = TRUE, scale = TRUE)
  }
  
  n <- length(y)
  if(Yscale) y <- scale(y, center=TRUE, scale=TRUE)
  
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
  
  outg <- optim(c(rep(1, ncol(X))), nlsep, method="L-BFGS-B", 
                lower=c(rep(0.1*sqrt(ncol(X)), ncol(X))), upper=c(rep(100000*sqrt(ncol(X)), ncol(X))), X=X, Y=y) 
  
  K <- covar.sep(X, d=outg$par, g=g)
  Ki <- solve(K+diag(g,n))
  tau2hat <- drop(t(y) %*% Ki %*% y / n)
  
  return(list(theta = outg$par, g=g, Ki=Ki, X = X, y = y, tau2hat=tau2hat, Xscale=Xscale, Yscale=Yscale))
}

pred.GP <- function(fit, xnew){
  
  xnew <- as.matrix(xnew)
  
  Xscale <- fit$Xscale
  Yscale <- fit$Yscale
  Ki <- fit$Ki
  theta <- fit$theta
  g <- fit$g
  X <- fit$X
  y <- fit$y
  tau2hat <- fit$tau2hat
  
  if(Xscale) xnew <- t((t(xnew)-attr(X,"scaled:center"))/attr(X,"scaled:scale"))
  
  KXX <- covar.sep(xnew, d=theta, g=g) 
  KX <- covar.sep(xnew, X, d=theta, g=0)
  
  if(Yscale) mup2 <- KX %*% Ki %*% (y * attr(y, "scaled:scale") + attr(y, "scaled:center")) else mup2 <- KX %*% Ki %*% y
  Sigmap2 <- pmax(0, diag(tau2hat*(KXX - KX %*% Ki %*% t(KX))))
  
  return(list(mu=mup2, sig2=Sigmap2))
}
