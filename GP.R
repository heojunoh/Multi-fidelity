library(plgp)
source("cor.sep.R")
source("matern.kernel.R")
eps <- sqrt(.Machine$double.eps) 

GP <- function(X, y, g=eps, 
               lower=0.001, upper=1000, 
               Xscale=TRUE, Yscale=TRUE, constant=FALSE){
  if(constant){
    
    if(is.null(dim(X))) X <- matrix(X, ncol = 1)
    
    if(Xscale){
      X <- scale(X, center = TRUE, scale = TRUE)
    }else{
      attr(X,"scaled:center") <- rep(0, ncol(X))
      attr(X,"scaled:scale") <- rep(1, ncol(X))
    }
    
    # darg way
    init <- rep(sort(distance(X))[sort(distance(X))!=0][0.1*length(sort(distance(X))[sort(distance(X))!=0])], ncol(X))
    # GPy way
    # init <- rep(median(distance(X)), ncol(X))
    # init <- rep(1, ncol(X))
    
    n <- length(y)
    
    nlsep <- function(par, X, Y) 
    {
      theta <- par # lengthscale
      K <- covar.sep(X, d=theta, g=g)
      Ki <- solve(K)
      ldetK <- determinant(K, logarithm=TRUE)$modulus
      
      one.vec <- matrix(1,ncol=1,nrow=n)
      mu.hat <- drop((t(one.vec)%*%Ki%*%Y)/(t(one.vec)%*%Ki%*%one.vec))
      
      tau2hat <- drop(t(Y-mu.hat)%*%Ki%*%(Y-mu.hat)/n)
      ll <- - (n/2)*log(tau2hat) - (1/2)*ldetK
      return(drop(-ll))
    }
    
    gradnlsep <- function(par, X, Y)
    {
      theta <- par
      K <- covar.sep(X, d=theta, g=g)
      Ki <- solve(K)
      
      one.vec <- matrix(1,ncol=1,nrow=n)
      mu.hat <- drop((t(one.vec)%*%Ki%*%Y)/(t(one.vec)%*%Ki%*%one.vec))
      
      KiY <- Ki %*% (Y-mu.hat)
      ## loop over theta components
      dlltheta <- rep(NA, length(theta))
      for(k in 1:length(dlltheta)){
        dotK <- K *distance(X[,k])/(theta[k]^2)
        dlltheta[k] <- (n/2) *t(KiY) %*% dotK %*% KiY / (t(Y) %*% KiY) - (1/2)*sum(diag(Ki %*% dotK))
      }

      return(-c(dlltheta))
    }
    
    outg <- optim(init, nlsep, gradnlsep, 
                  method="L-BFGS-B", lower=lower, upper=upper, X=X, Y=y) 
    
    K <- covar.sep(X, d=outg$par, g=g)
    Ki <- solve(K)
    one.vec <- matrix(1,ncol=1,nrow=n)
    mu.hat <- drop((t(one.vec)%*%Ki%*%y)/(t(one.vec)%*%Ki%*%one.vec))
    tau2hat <- drop(t(y-mu.hat) %*% Ki %*% (y-mu.hat) / nrow(X))
    
    return(list(theta = outg$par, g=g, Ki=Ki, mu.hat=mu.hat, X = X, y = y, tau2hat=tau2hat, Xscale=Xscale, constant=constant))
  }else{
    
    if(is.null(dim(X))) X <- matrix(X, ncol = 1)
    
    if(Xscale){
      X <- scale(X, center = TRUE, scale = TRUE)
    }else{
      attr(X,"scaled:center") <- rep(0, ncol(X))
      attr(X,"scaled:scale") <- rep(1, ncol(X))
    }
    
    # darg way
    init <- rep(sort(distance(X))[sort(distance(X))!=0][0.1*length(sort(distance(X))[sort(distance(X))!=0])], ncol(X))
    # GPy way
    # init <- rep(median(distance(X)), ncol(X))
    
    n <- length(y)
    if(Yscale) y <- scale(y, center=TRUE, scale=FALSE) # If use mean, don't scale
    
    nlsep <- function(par, X, Y) 
    {
      theta <- par # lengthscale
      K <- covar.sep(X, d=theta, g=g)
      Ki <- solve(K)
      ldetK <- determinant(K, logarithm=TRUE)$modulus
      ll <- - (n/2)*log(t(Y) %*% Ki %*% Y) - (1/2)*ldetK
      return(drop(-ll))
    }
    
    # gradnlsep <- function(par, X, Y)
    # {
    #   theta <- par
    #   n <- length(Y) 
    #   K <- covar.sep(X, d=theta, g=g) 
    #   Ki <- solve(K) 
    #   KiY <- Ki %*% Y 
    #   
    #   ## loop over theta components 
    #   dlltheta <- rep(NA, length(theta)) 
    #   for(k in 1:length(dlltheta)){ 
    #     dotK <- K *distance(X[,k])/(theta[k]^2) 
    #     dlltheta[k] <- (n/2) *t(KiY) %*% dotK %*% KiY / (t(Y) %*% KiY) - (1/2)*sum(diag(Ki %*% dotK)) 
    #   } 
    #   
    #   return(-c(dlltheta)) 
    # }
    
    outg <- optim(init, nlsep, #gradnlsep, 
                  method="L-BFGS-B", lower=lower, upper=upper, X=X, Y=y) 

    K <- covar.sep(X, d=outg$par, g=g)
    Ki <- solve(K)
    tau2hat <- drop(t(y) %*% Ki %*% y / n)
    mu.hat <- 0
    
    return(list(theta = outg$par, g=g, Ki=Ki, mu.hat=mu.hat, X = X, y = y, tau2hat=tau2hat, Xscale=Xscale, Yscale=Yscale, constant=constant))
  }
}

matGP <- function(X, y, nu=2.5, g=eps, 
                  lower=rep(0.1, ncol(X)), upper=rep(100,ncol(X)), 
                  Xscale=TRUE, Yscale=TRUE, constant=FALSE, indicator=NULL){
  if(constant){
    
    if(is.null(dim(X))) X <- matrix(X, ncol = 1)
    
    if(Xscale){
      X <- scale(X, center = TRUE, scale = TRUE)
    }else{
      attr(X,"scaled:center") <- rep(0, ncol(X))
      attr(X,"scaled:scale") <- rep(1, ncol(X))
    }
    
    init <- rep(1*sqrt(ncol(X)), ncol(X))
    n <- length(y)
    
    nlsep <- function(par, X, Y) 
    {
      theta <- par # lengthscale
      K <- cor.sep(X, theta=theta, nu=nu)
      # R <- sqrt(distance(t(t(X)/theta)))
      # K <- matern.kernel(R, nu=nu)
      Ki <- solve(K+diag(g,n))
      ldetK <- determinant(K, logarithm=TRUE)$modulus
      
      one.vec <- matrix(1,ncol=1,nrow=n)
      mu.hat <- drop((t(one.vec)%*%Ki%*%Y)/(t(one.vec)%*%Ki%*%one.vec))
      
      tau2hat <- drop(t(Y-mu.hat)%*%Ki%*%(Y-mu.hat)/n)
      ll <- - (n/2)*log(tau2hat) - (1/2)*ldetK
      return(-ll)
    }
    
    outg <- optim(init, nlsep,  
                  method="L-BFGS-B", lower=lower, upper=upper, X=X, Y=y) 
    
    theta <- outg$par
    
    R <- sqrt(distance(t(t(X)/theta)))
    K <- matern.kernel(R, nu=nu)
    Ki <- solve(K+diag(g,n))
    one.vec <- matrix(1,ncol=1,nrow=n)
    mu.hat <- drop((t(one.vec)%*%Ki%*%y)/(t(one.vec)%*%Ki%*%one.vec))
    tau2hat <- drop(t(y-mu.hat) %*% Ki %*% (y-mu.hat) / nrow(X))
    
    return(list(theta = theta, nu=nu, g=g, Ki=Ki, mu.hat=mu.hat, X = X, y = y, tau2hat=tau2hat, Xscale=Xscale, constant=constant))
  }else{
    
    if(is.null(dim(X))) X <- matrix(X, ncol = 1)
    
    if(Xscale){
      X <- scale(X, center = TRUE, scale = TRUE)
    }else{
      attr(X,"scaled:center") <- rep(0, ncol(X))
      attr(X,"scaled:scale") <- rep(1, ncol(X))
    }
    
    init <- rep(1*sqrt(ncol(X)), ncol(X))
    n <- length(y)
    if(Yscale) y <- scale(y, center=TRUE, scale=FALSE) # If use mean, don't scale
    
    nlsep <- function(par, X, Y) 
    {
      theta <- par # lengthscale
      K <- cor.sep(X, theta=theta, nu=nu)
      Ki <- solve(K+diag(g,n))
      ldetK <- determinant(K, logarithm=TRUE)$modulus
      tau2hat <- drop(t(Y)%*%Ki%*%(Y)/n)
      ll <- - (n/2)*log(tau2hat) - (1/2)*ldetK
      return(drop(-ll))
    }
    
    outg <- optim(init, nlsep, #gradnlsep, 
                  method="L-BFGS-B", lower=lower, upper=upper, X=X, Y=y) 
    
    K <- cor.sep(X, theta=outg$par, nu=nu)
    Ki <- solve(K+diag(g,n))
    tau2hat <- drop(t(y) %*% Ki %*% (y) / n)
    mu.hat <- 0
    
    return(list(theta = outg$par, nu=nu, g=g, Ki=Ki, mu.hat=mu.hat, X = X, y = y, tau2hat=tau2hat, Xscale=Xscale, Yscale=Yscale, constant=constant))
  }
}



pred.GP <- function(fit, xnew){
  constant <- fit$constant
  
  if(constant){
    xnew <- as.matrix(xnew)
    
    Xscale <- fit$Xscale
    Ki <- fit$Ki
    theta <- fit$theta
    g <- fit$g
    X <- fit$X
    y <- fit$y
    tau2hat <- fit$tau2hat
    mu.hat <- fit$mu.hat
    
    if(Xscale) xnew <- t((t(xnew)-attr(X,"scaled:center"))/attr(X,"scaled:scale"))
    
    KXX <- covar.sep(xnew, d=theta, g=g) 
    KX <- covar.sep(xnew, X, d=theta, g=0)
    
    mup2 <- mu.hat + KX %*% Ki %*% (y - mu.hat)
    Sigmap2 <- pmax(0, diag(tau2hat*(KXX - KX %*% Ki %*% t(KX))))
    
    return(list(mu=mup2, sig2=Sigmap2))
  }else{
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
    
    if(Yscale) mup2 <- KX %*% Ki %*% (y + attr(y, "scaled:center")) else mup2 <- KX %*% Ki %*% y
    Sigmap2 <- pmax(0, diag(tau2hat*(KXX - KX %*% Ki %*% t(KX))))
    
    return(list(mu=mup2, sig2=Sigmap2))
  }
}

pred.matGP <- function(fit, xnew){
  constant <- fit$constant
  
  if(constant){
    xnew <- as.matrix(xnew)
    
    Xscale <- fit$Xscale
    Ki <- fit$Ki
    theta <- fit$theta
    nu <- fit$nu
    g <- fit$g
    X <- fit$X
    y <- fit$y
    tau2hat <- fit$tau2hat
    mu.hat <- fit$mu.hat
    
    if(Xscale) xnew <- t((t(xnew)-attr(X,"scaled:center"))/attr(X,"scaled:scale"))
    
    KXX <- cor.sep(xnew, theta=theta, nu=nu)
    KX <- t(cor.sep(X, xnew, theta=theta, nu=nu))
    # RXX <- sqrt(distance(t(t(xnew)/theta)))
    # RX <- sqrt(distance(t(t(xnew)/theta), t(t(X)/theta)))
    # KXX <- matern.kernel(RXX, nu=nu)
    # KX <- matern.kernel(RX, nu=nu)
    
    mup2 <- mu.hat + KX %*% Ki %*% (y - mu.hat)
    Sigmap2 <- pmax(0, diag(tau2hat*(KXX + diag(g,nrow(xnew)) - KX %*% Ki %*% t(KX))))
    
    return(list(mu=mup2, sig2=Sigmap2))
  }else{
    xnew <- as.matrix(xnew)
    
    Xscale <- fit$Xscale
    Yscale <- fit$Yscale
    Ki <- fit$Ki
    theta <- fit$theta
    nu <- fit$nu
    g <- fit$g
    X <- fit$X
    y <- fit$y
    tau2hat <- fit$tau2hat
    
    if(Xscale) xnew <- t((t(xnew)-attr(X,"scaled:center"))/attr(X,"scaled:scale"))
    
    KXX <- cor.sep(xnew, theta=theta, nu=nu)
    KX <- t(cor.sep(X, xnew, theta=theta, nu=nu))
    
    if(Yscale) mup2 <- KX %*% Ki %*% (y + attr(y, "scaled:center")) else mup2 <- KX %*% Ki %*% y
    Sigmap2 <- pmax(0, diag(tau2hat*(KXX + diag(g,nrow(xnew)) - KX %*% Ki %*% t(KX))))
    
    return(list(mu=mup2, sig2=Sigmap2))
  }
}

