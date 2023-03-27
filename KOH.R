library(plgp)
eps <- sqrt(.Machine$double.eps) 

KOHGP <- function(X, y, g=eps, center=TRUE){
  
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
  
  outg <- optim(c(rep(1, ncol(X))), nlsep, method="L-BFGS-B", 
                lower=c(rep(0.1*sqrt(ncol(X)), ncol(X))), upper=c(rep(100000*sqrt(ncol(X)), ncol(X))), X=X, Y=y) 
  
  K <- covar.sep(X, d=outg$par, g=g)
  Ki <- solve(K+diag(g,n))
  tau2hat <- drop(t(y) %*% Ki %*% y / n)
  
  return(list(theta = outg$par, g=g, Ki=Ki, X = X, y = y, tau2hat=tau2hat))
}

KOH <- function(X, Y2, Y1, g=eps, center=TRUE){
  
  n <- length(Y2)
  # if(center) y <- scale(y, scale=FALSE)
  
  nlsep2 <- function(par, X, Y2, Y1) 
  {
    n <- length(Y2)
    theta <- par # lengthscale and rho
    Y <- Y2-theta[ncol(X)+1]*Y1 # y2-rho*y1
    K <- covar.sep(X, d=theta[1:ncol(X)], g=g)
    Ki <- solve(K+diag(g,n))
    ldetK <- determinant(K, logarithm=TRUE)$modulus
    ll <- - (n/2)*log(t(Y) %*% Ki %*% Y) - (1/2)*ldetK
    return(drop(-ll))
  }
  
  outg <- optim(c(rep(1, ncol(X)), 0.1), # For lengthscale and rho #
                nlsep2, method="L-BFGS-B", 
                lower=c(rep(0.1*sqrt(ncol(X)), ncol(X)+1)), 
                upper=c(rep(100000*sqrt(ncol(X)), ncol(X)+1)), 
                X=X, Y2=Y2, Y1=Y1) 
  
  K <- covar.sep(X, d=outg$par[1:ncol(X)], g=g)
  Ki <- solve(K+diag(g,n))
  rho <- outg$par[ncol(X)+1]
  y <- Y2-rho*Y1
  tau2hat <- drop(t(y) %*% Ki %*% y / n)
  
  return(list(theta = outg$par[1:ncol(X)], rho = rho, g=g, Ki=Ki, X = X, y = y, tau2hat=tau2hat))
}

fit.KOH <- function(X1, X2, X3, Y1, Y2, Y3, g=eps){ # need to change function for another example
  
  ### KOH method ###
  Y2d3 <- fl(X3, l=3)
  Y1d2 <- fl(X2, l=1)
  
  ### estimating first order ###
  fit.KOHGP1 <- KOHGP(X1, Y1)
  b1 <- 1/fit.KOHGP1$theta
  sig2_1 <- fit.KOHGP1$tau2hat
  
  ### estimating second order ###
  # KOH(X2, Y2, Y1d2)
  rho1 <- KOH(X2, Y2, Y1d2)$rho
  b2 <- 1/KOH(X2, Y2, Y1d2)$theta
  sig2_2 <- KOH(X2, Y2, Y1d2)$tau2hat
  
  ### estimating third order ###
  # KOH(X3, Y3, Y2d3)
  rho2 <- KOH(X3, Y3, Y2d3)$rho
  b3 <- 1/KOH(X3, Y3, Y2d3)$theta
  sig2_3 <- KOH(X3, Y3, Y2d3)$tau2hat
  
  return(list(b=c(b1, b2, b3), rho=c(rho1, rho2), tau2hat=c(sig2_1, sig2_2, sig2_3), g=g, X1=X1, X2=X2, X3=X3, Y1=Y1, Y2=Y2, Y3=Y3))
}

pred.KOH <- function(fit, x){ # need to change function for another example
  
  X1 <- fit$X1
  X2 <- fit$X2
  X3 <- fit$X3
  Y1 <- fit$Y1
  Y2 <- fit$Y2
  Y3 <- fit$Y3
  
  b <- fit$b
  rho <- fit$rho
  tau2hat <- fit$tau2hat
  g <- fit$g
  
  
  ### prediction of 2nd order KOH ###
  tx2 <- cbind(rho[1]*rho[2]*tau2hat[1]*covar.sep(x, X1, d=1/b[1], g=g), 
               rho[1]^2*rho[2]*tau2hat[1]*covar.sep(x, X2, d=1/b[1], g=g) + rho[2]*tau2hat[2]*covar.sep(x, X2, d=1/b[2], g=g),
               rho[1]^2*rho[2]^2*tau2hat[1]*covar.sep(x, X3, d=1/b[1], g=g) + rho[2]^2*tau2hat[2]*covar.sep(x, X3, d=1/b[2], g=g) + tau2hat[3]*covar.sep(x, X3, d=1/b[3], g=g))
  
  V1 <- tau2hat[1]*covar.sep(X1, d=1/b[1], g=g)
  V12 <- rho[1]*tau2hat[1]*covar.sep(X1, X2, d=1/b[1], g=0)
  V13 <- rho[1]*rho[2]*tau2hat[1]*covar.sep(X1, X3, d=1/b[1], g=0)
  V2 <- rho[1]^2*tau2hat[1]*covar.sep(X2, d=1/b[1], g=g) + tau2hat[2]*covar.sep(X2, d=1/b[2], g=g)
  V23 <- rho[1]^2*rho[2]*tau2hat[1]*covar.sep(X2, X3, d=1/b[1], g=0) + rho[2]*tau2hat[2]*covar.sep(X2, X3, d=1/b[2], g=0)
  V3 <- rho[1]^2*rho[2]^2*tau2hat[1]*covar.sep(X3, d=1/b[1], g=g) + rho[2]^2*tau2hat[2]*covar.sep(X3, d=1/b[2], g=g) + tau2hat[3]*covar.sep(X3, d=1/b[3], g=g)
  
  V_3 <- rbind(cbind(V1, V12, V13), cbind(t(V12), V2, V23), cbind(t(V13), t(V23), V3))
  
  mx2 <- tx2 %*% solve(V_3) %*% c(Y1, Y2, Y3)
  
  ### posterior variance ###
  koh.var2 <- pmax(0, diag(tau2hat[3]*covar.sep(matrix(x), d=1/b[3], g=g) + tau2hat[2]*rho[2]^2*covar.sep(matrix(x), d=1/b[2], g=g) + tau2hat[1]*rho[1]^2*rho[2]^2*covar.sep(matrix(x), d=1/b[3], g=g) - tx2 %*% solve(V_3+diag(g, nrow(V_3)))%*%t(tx2)))
  
  return(list(mu=mx2, sig2=koh.var2))
}


IMSPEKOH1 <- function(x, newx, fit, mc.sample=10000){
  ### This is updating KOH when the one data point is added ###
  # x; usually grid
  # newx; new point which will be added, should be n=1
  # fit; fitted KOH model
  X1 <- fit$X1
  X2 <- fit$X2
  X3 <- fit$X3
  Y1 <- fit$Y1
  Y2 <- fit$Y2
  Y3 <- fit$Y3
  
  b <- fit$b
  rho <- fit$rho
  tau2hat <- fit$tau2hat
  g <- fit$g
  
  ### Generate MC output
  m1 <- covar.sep(newx, X1, d=1/b[1], g=0)  %*% solve(covar.sep(X1, d=1/b[1], g=g)) %*% Y1
  v1 <- pmax(0, diag(tau2hat[1]*(covar.sep(matrix(newx), d=1/b[1], g=g) -
                                   covar.sep(newx, X1, d=1/b[1], g=0) %*% solve(covar.sep(X1, d=1/b[1], g=g)) %*% t(covar.sep(newx, X1, d=1/b[1], g=0)))))
  
  x1.sample <- rnorm(mc.sample, mean=m1, sd=sqrt(v1))
  mu.cand1 <- mean(x1.sample) # f1(newx)
  
  X1 <- rbind(fit$X1, newx)
  Y1 <- rbind(fit$Y1, mu.cand1)
  
  ### update sig2
  tx2 <- cbind(rho[1]*rho[2]*tau2hat[1]*covar.sep(x, X1, d=1/b[1], g=g), 
               rho[1]^2*rho[2]*tau2hat[1]*covar.sep(x, X2, d=1/b[1], g=g) + rho[2]*tau2hat[2]*covar.sep(x, X2, d=1/b[2], g=g),
               rho[1]^2*rho[2]^2*tau2hat[1]*covar.sep(x, X3, d=1/b[1], g=g) + rho[2]^2*tau2hat[2]*covar.sep(x, X3, d=1/b[2], g=g) + tau2hat[3]*covar.sep(x, X3, d=1/b[3], g=g))
  
  V1 <- tau2hat[1]*covar.sep(X1, d=1/b[1], g=g)
  V12 <- rho[1]*tau2hat[1]*covar.sep(X1, X2, d=1/b[1], g=0)
  V13 <- rho[1]*rho[2]*tau2hat[1]*covar.sep(X1, X3, d=1/b[1], g=0)
  V2 <- rho[1]^2*tau2hat[1]*covar.sep(X2, d=1/b[1], g=g) + tau2hat[2]*covar.sep(X2, d=1/b[2], g=g)
  V23 <- rho[1]^2*rho[2]*tau2hat[1]*covar.sep(X2, X3, d=1/b[1], g=0) + rho[2]*tau2hat[2]*covar.sep(X2, X3, d=1/b[2], g=0)
  V3 <- rho[1]^2*rho[2]^2*tau2hat[1]*covar.sep(X3, d=1/b[1], g=g) + rho[2]^2*tau2hat[2]*covar.sep(X3, d=1/b[2], g=g) + tau2hat[3]*covar.sep(X3, d=1/b[3], g=g)
  
  V_3 <- rbind(cbind(V1, V12, V13), cbind(t(V12), V2, V23), cbind(t(V13), t(V23), V3))
  
  mx2 <- tx2 %*% solve(V_3) %*% c(Y1, Y2, Y3)
  
  koh.var2 <- pmax(0, diag(tau2hat[3]*covar.sep(matrix(x), d=1/b[3], g=g) + tau2hat[2]*rho[2]^2*covar.sep(matrix(x), d=1/b[2], g=g) + tau2hat[1]*rho[1]^2*rho[2]^2*covar.sep(matrix(x), d=1/b[3], g=g) - tx2 %*% solve(V_3+diag(g, nrow(V_3)))%*%t(tx2)))
  
  return(IMSPE=mean(koh.var2))
}


IMSPEKOH2 <- function(x, newx, fit, mc.sample=10000){
  ### This is updating KOH when the one data point is added ###
  # x; usually grid
  # newx; new point which will be added, should be n=1
  # fit; fitted KOH model
  X1 <- fit$X1
  X2 <- fit$X2
  X3 <- fit$X3
  Y1 <- fit$Y1
  Y2 <- fit$Y2
  Y3 <- fit$Y3
  
  b <- fit$b
  rho <- fit$rho
  tau2hat <- fit$tau2hat
  g <- fit$g
  
  ### Generate MC output
  m1 <- covar.sep(newx, X1, d=1/b[1], g=0)  %*% solve(covar.sep(X1, d=1/b[1], g=g)) %*% Y1
  var1 <- pmax(0, diag(tau2hat[1]*(covar.sep(matrix(newx), d=1/b[1], g=g) -
                                     covar.sep(newx, X1, d=1/b[1], g=0) %*% solve(covar.sep(X1, d=1/b[1], g=g)) %*% t(covar.sep(newx, X1, d=1/b[1], g=0)))))
  
  x1.sample <- rnorm(mc.sample, mean=m1, sd=sqrt(var1))
  mu.cand1 <- mean(x1.sample) # f1(newx)
  
  X1 <- rbind(fit$X1, newx)
  Y1 <- rbind(fit$Y1, mu.cand1)  
  
  ### 2nd order KOH ###
  tx1 <- cbind(rho[1]*tau2hat[1]*covar.sep(newx, X1, d=1/b[1], g=g), 
               rho[1]^2*tau2hat[1]*covar.sep(newx, X2, d=1/b[1], g=g) + tau2hat[2]*covar.sep(newx, X2, d=1/b[2], g=g))
  
  V1 <- tau2hat[1]*covar.sep(X1, d=1/b[1], g=g)
  V12 <- rho[1]*tau2hat[1]*covar.sep(X1, X2, d=1/b[1], g=0)
  V2 <- rho[1]^2*tau2hat[1]*covar.sep(X2, d=1/b[1], g=g) + tau2hat[2]*covar.sep(X2, d=1/b[2], g=g)
  
  V_2 <- rbind(cbind(V1, V12), cbind(t(V12), V2))
  
  m2 <- tx1 %*% solve(V_2) %*% c(Y1, Y2)
  var2 <- pmax(0, diag(tau2hat[2]*covar.sep(matrix(newx), d=1/b[2], g=g) + tau2hat[1]*rho[1]^2*covar.sep(matrix(newx), d=1/b[1], g=g) - tx1 %*% solve(V_2+diag(g, nrow(V_2)))%*%t(tx1)))
  
  x2.sample <- rnorm(mc.sample, mean=m2, sd=sqrt(var2))
  mu.cand2 <- mean(x2.sample) # f2(newx)
  
  X2 <- rbind(fit$X2, newx)
  Y2 <- rbind(fit$Y2, mu.cand2)
  
  ### update sig2
  tx2 <- cbind(rho[1]*rho[2]*tau2hat[1]*covar.sep(x, X1, d=1/b[1], g=g), 
               rho[1]^2*rho[2]*tau2hat[1]*covar.sep(x, X2, d=1/b[1], g=g) + rho[2]*tau2hat[2]*covar.sep(x, X2, d=1/b[2], g=g),
               rho[1]^2*rho[2]^2*tau2hat[1]*covar.sep(x, X3, d=1/b[1], g=g) + rho[2]^2*tau2hat[2]*covar.sep(x, X3, d=1/b[2], g=g) + tau2hat[3]*covar.sep(x, X3, d=1/b[3], g=g))
  
  V1 <- tau2hat[1]*covar.sep(X1, d=1/b[1], g=g)
  V12 <- rho[1]*tau2hat[1]*covar.sep(X1, X2, d=1/b[1], g=0)
  V13 <- rho[1]*rho[2]*tau2hat[1]*covar.sep(X1, X3, d=1/b[1], g=0)
  V2 <- rho[1]^2*tau2hat[1]*covar.sep(X2, d=1/b[1], g=g) + tau2hat[2]*covar.sep(X2, d=1/b[2], g=g)
  V23 <- rho[1]^2*rho[2]*tau2hat[1]*covar.sep(X2, X3, d=1/b[1], g=0) + rho[2]*tau2hat[2]*covar.sep(X2, X3, d=1/b[2], g=0)
  V3 <- rho[1]^2*rho[2]^2*tau2hat[1]*covar.sep(X3, d=1/b[1], g=g) + rho[2]^2*tau2hat[2]*covar.sep(X3, d=1/b[2], g=g) + tau2hat[3]*covar.sep(X3, d=1/b[3], g=g)
  
  V_3 <- rbind(cbind(V1, V12, V13), cbind(t(V12), V2, V23), cbind(t(V13), t(V23), V3))
  
  mx2 <- tx2 %*% solve(V_3) %*% c(Y1, Y2, Y3)
  
  koh.var2 <- pmax(0, diag(tau2hat[3]*covar.sep(matrix(x), d=1/b[3], g=g) + tau2hat[2]*rho[2]^2*covar.sep(matrix(x), d=1/b[2], g=g) + tau2hat[1]*rho[1]^2*rho[2]^2*covar.sep(matrix(x), d=1/b[3], g=g) - tx2 %*% solve(V_3+diag(g, nrow(V_3)))%*%t(tx2)))
  
  return(IMSPE=mean(koh.var2))
}


IMSPEKOH3 <- function(x, newx, fit, mc.sample=10000){
  ### This is updating KOH when the one data point is added ###
  # x; usually grid
  # newx; new point which will be added, should be n=1
  # fit; fitted KOH model
  X1 <- fit$X1
  X2 <- fit$X2
  X3 <- fit$X3
  Y1 <- fit$Y1
  Y2 <- fit$Y2
  Y3 <- fit$Y3
  
  b <- fit$b
  rho <- fit$rho
  tau2hat <- fit$tau2hat
  g <- fit$g
  
  ### Generate MC output
  m1 <- covar.sep(newx, X1, d=1/b[1], g=0)  %*% solve(covar.sep(X1, d=1/b[1], g=g)) %*% Y1
  var1 <- pmax(0, diag(tau2hat[1]*(covar.sep(matrix(newx), d=1/b[1], g=g) -
                                     covar.sep(newx, X1, d=1/b[1], g=0) %*% solve(covar.sep(X1, d=1/b[1], g=g)) %*% t(covar.sep(newx, X1, d=1/b[1], g=0)))))
  
  x1.sample <- rnorm(mc.sample, mean=m1, sd=sqrt(var1))
  mu.cand1 <- mean(x1.sample) # f1(newx)
  
  X1 <- rbind(fit$X1, newx)
  Y1 <- rbind(fit$Y1, mu.cand1)  
  
  ### 2nd order KOH ###
  tx1 <- cbind(rho[1]*tau2hat[1]*covar.sep(newx, X1, d=1/b[1], g=g), 
               rho[1]^2*tau2hat[1]*covar.sep(newx, X2, d=1/b[1], g=g) + tau2hat[2]*covar.sep(newx, X2, d=1/b[2], g=g))
  
  V1 <- tau2hat[1]*covar.sep(X1, d=1/b[1], g=g)
  V12 <- rho[1]*tau2hat[1]*covar.sep(X1, X2, d=1/b[1], g=0)
  V2 <- rho[1]^2*tau2hat[1]*covar.sep(X2, d=1/b[1], g=g) + tau2hat[2]*covar.sep(X2, d=1/b[2], g=g)
  
  V_2 <- rbind(cbind(V1, V12), cbind(t(V12), V2))
  
  m2 <- tx1 %*% solve(V_2) %*% c(Y1, Y2)
  var2 <- pmax(0, diag(tau2hat[2]*covar.sep(matrix(newx), d=1/b[2], g=g) + tau2hat[1]*rho[1]^2*covar.sep(matrix(newx), d=1/b[1], g=g) - tx1 %*% solve(V_2+diag(g, nrow(V_2)))%*%t(tx1)))
  
  x2.sample <- rnorm(mc.sample, mean=m2, sd=sqrt(var2))
  mu.cand2 <- mean(x2.sample) # f2(newx)
  
  X2 <- rbind(fit$X2, newx)
  Y2 <- rbind(fit$Y2, mu.cand2)
  
  ### 3rd order KOH ###
  tx2 <- cbind(rho[1]*rho[2]*tau2hat[1]*covar.sep(newx, X1, d=1/b[1], g=g), 
               rho[1]^2*rho[2]*tau2hat[1]*covar.sep(newx, X2, d=1/b[1], g=g) + rho[2]*tau2hat[2]*covar.sep(newx, X2, d=1/b[2], g=g),
               rho[1]^2*rho[2]^2*tau2hat[1]*covar.sep(newx, X3, d=1/b[1], g=g) + rho[2]^2*tau2hat[2]*covar.sep(newx, X3, d=1/b[2], g=g) + tau2hat[3]*covar.sep(newx, X3, d=1/b[3], g=g))
  
  V1 <- tau2hat[1]*covar.sep(X1, d=1/b[1], g=g)
  V12 <- rho[1]*tau2hat[1]*covar.sep(X1, X2, d=1/b[1], g=0)
  V13 <- rho[1]*rho[2]*tau2hat[1]*covar.sep(X1, X3, d=1/b[1], g=0)
  V2 <- rho[1]^2*tau2hat[1]*covar.sep(X2, d=1/b[1], g=g) + tau2hat[2]*covar.sep(X2, d=1/b[2], g=g)
  V23 <- rho[1]^2*rho[2]*tau2hat[1]*covar.sep(X2, X3, d=1/b[1], g=0) + rho[2]*tau2hat[2]*covar.sep(X2, X3, d=1/b[2], g=0)
  V3 <- rho[1]^2*rho[2]^2*tau2hat[1]*covar.sep(X3, d=1/b[1], g=g) + rho[2]^2*tau2hat[2]*covar.sep(X3, d=1/b[2], g=g) + tau2hat[3]*covar.sep(X3, d=1/b[3], g=g)
  
  V_3 <- rbind(cbind(V1, V12, V13), cbind(t(V12), V2, V23), cbind(t(V13), t(V23), V3))
  
  m3 <- tx2 %*% solve(V_3) %*% c(Y1, Y2, Y3)
  
  var3 <- pmax(0, diag(tau2hat[3]*covar.sep(matrix(newx), d=1/b[3], g=g) + tau2hat[2]*rho[2]^2*covar.sep(matrix(newx), d=1/b[2], g=g) + tau2hat[1]*rho[1]^2*rho[2]^2*covar.sep(matrix(newx), d=1/b[3], g=g) - tx2 %*% solve(V_3+diag(g, nrow(V_3)))%*%t(tx2)))
  
  x3.sample <- rnorm(mc.sample, mean=m3, sd=sqrt(var3))
  mu.cand3 <- mean(x3.sample) # f2(newx)
  
  X3 <- rbind(fit$X3, newx)
  Y3 <- rbind(fit$Y3, mu.cand3)
  
  ### update sig2
  tx2 <- cbind(rho[1]*rho[2]*tau2hat[1]*covar.sep(x, X1, d=1/b[1], g=g), 
               rho[1]^2*rho[2]*tau2hat[1]*covar.sep(x, X2, d=1/b[1], g=g) + rho[2]*tau2hat[2]*covar.sep(x, X2, d=1/b[2], g=g),
               rho[1]^2*rho[2]^2*tau2hat[1]*covar.sep(x, X3, d=1/b[1], g=g) + rho[2]^2*tau2hat[2]*covar.sep(x, X3, d=1/b[2], g=g) + tau2hat[3]*covar.sep(x, X3, d=1/b[3], g=g))
  
  V1 <- tau2hat[1]*covar.sep(X1, d=1/b[1], g=g)
  V12 <- rho[1]*tau2hat[1]*covar.sep(X1, X2, d=1/b[1], g=0)
  V13 <- rho[1]*rho[2]*tau2hat[1]*covar.sep(X1, X3, d=1/b[1], g=0)
  V2 <- rho[1]^2*tau2hat[1]*covar.sep(X2, d=1/b[1], g=g) + tau2hat[2]*covar.sep(X2, d=1/b[2], g=g)
  V23 <- rho[1]^2*rho[2]*tau2hat[1]*covar.sep(X2, X3, d=1/b[1], g=0) + rho[2]*tau2hat[2]*covar.sep(X2, X3, d=1/b[2], g=0)
  V3 <- rho[1]^2*rho[2]^2*tau2hat[1]*covar.sep(X3, d=1/b[1], g=g) + rho[2]^2*tau2hat[2]*covar.sep(X3, d=1/b[2], g=g) + tau2hat[3]*covar.sep(X3, d=1/b[3], g=g)
  
  V_3 <- rbind(cbind(V1, V12, V13), cbind(t(V12), V2, V23), cbind(t(V13), t(V23), V3))
  
  mx2 <- tx2 %*% solve(V_3) %*% c(Y1, Y2, Y3)
  
  koh.var2 <- pmax(0, diag(tau2hat[3]*covar.sep(matrix(x), d=1/b[3], g=g) + tau2hat[2]*rho[2]^2*covar.sep(matrix(x), d=1/b[2], g=g) + tau2hat[1]*rho[1]^2*rho[2]^2*covar.sep(matrix(x), d=1/b[3], g=g) - tx2 %*% solve(V_3+diag(g, nrow(V_3)))%*%t(tx2)))
  
  return(IMSPE=mean(koh.var2))
}


IMSPEKOH1select <- function(x, newx, fit, mc.sample=10000){
  ### This is updating KOH when the one data point is added ###
  # x; usually grid
  # newx; new point which will be added, should be n=1
  # fit; fitted KOH model
  X1 <- fit$X1
  X2 <- fit$X2
  X3 <- fit$X3
  Y1 <- fit$Y1
  Y2 <- fit$Y2
  Y3 <- fit$Y3
  
  b <- fit$b
  rho <- fit$rho
  tau2hat <- fit$tau2hat
  g <- fit$g
  
  ### Generate output
  y1.select <- fl(newx, l=1)
  
  X1 <- rbind(fit$X1, newx)
  Y1 <- rbind(fit$Y1, y1.select)
  
  ### update sig2
  tx2 <- cbind(rho[1]*rho[2]*tau2hat[1]*covar.sep(x, X1, d=1/b[1], g=g), 
               rho[1]^2*rho[2]*tau2hat[1]*covar.sep(x, X2, d=1/b[1], g=g) + rho[2]*tau2hat[2]*covar.sep(x, X2, d=1/b[2], g=g),
               rho[1]^2*rho[2]^2*tau2hat[1]*covar.sep(x, X3, d=1/b[1], g=g) + rho[2]^2*tau2hat[2]*covar.sep(x, X3, d=1/b[2], g=g) + tau2hat[3]*covar.sep(x, X3, d=1/b[3], g=g))
  
  V1 <- tau2hat[1]*covar.sep(X1, d=1/b[1], g=g)
  V12 <- rho[1]*tau2hat[1]*covar.sep(X1, X2, d=1/b[1], g=0)
  V13 <- rho[1]*rho[2]*tau2hat[1]*covar.sep(X1, X3, d=1/b[1], g=0)
  V2 <- rho[1]^2*tau2hat[1]*covar.sep(X2, d=1/b[1], g=g) + tau2hat[2]*covar.sep(X2, d=1/b[2], g=g)
  V23 <- rho[1]^2*rho[2]*tau2hat[1]*covar.sep(X2, X3, d=1/b[1], g=0) + rho[2]*tau2hat[2]*covar.sep(X2, X3, d=1/b[2], g=0)
  V3 <- rho[1]^2*rho[2]^2*tau2hat[1]*covar.sep(X3, d=1/b[1], g=g) + rho[2]^2*tau2hat[2]*covar.sep(X3, d=1/b[2], g=g) + tau2hat[3]*covar.sep(X3, d=1/b[3], g=g)
  
  V_3 <- rbind(cbind(V1, V12, V13), cbind(t(V12), V2, V23), cbind(t(V13), t(V23), V3))
  
  mx2 <- tx2 %*% solve(V_3) %*% c(Y1, Y2, Y3)
  
  koh.var2 <- pmax(0, diag(tau2hat[3]*covar.sep(matrix(x), d=1/b[3], g=g) + tau2hat[2]*rho[2]^2*covar.sep(matrix(x), d=1/b[2], g=g) + tau2hat[1]*rho[1]^2*rho[2]^2*covar.sep(matrix(x), d=1/b[3], g=g) - tx2 %*% solve(V_3+diag(g, nrow(V_3)))%*%t(tx2)))
  
  fit$X1 <- X1
  fit$Y1 <- Y1
  
  return(list(IMSPE=mean(koh.var2), fit=fit))
}


IMSPEKOH2select <- function(x, newx, fit, mc.sample=10000){
  ### This is updating KOH when the one data point is added ###
  # x; usually grid
  # newx; new point which will be added, should be n=1
  # fit; fitted KOH model
  X1 <- fit$X1
  X2 <- fit$X2
  X3 <- fit$X3
  Y1 <- fit$Y1
  Y2 <- fit$Y2
  Y3 <- fit$Y3
  
  b <- fit$b
  rho <- fit$rho
  tau2hat <- fit$tau2hat
  g <- fit$g
  
  ### Generate output
  y1.select <- fl(newx, l=1)
  
  X1 <- rbind(fit$X1, newx)
  Y1 <- rbind(fit$Y1, y1.select)  
  
  ### 2nd order KOH ###
  y2.select <- fl(newx, l=3)
  
  X2 <- rbind(fit$X2, newx)
  Y2 <- rbind(fit$Y2, y2.select)
  
  ### update sig2
  tx2 <- cbind(rho[1]*rho[2]*tau2hat[1]*covar.sep(x, X1, d=1/b[1], g=g), 
               rho[1]^2*rho[2]*tau2hat[1]*covar.sep(x, X2, d=1/b[1], g=g) + rho[2]*tau2hat[2]*covar.sep(x, X2, d=1/b[2], g=g),
               rho[1]^2*rho[2]^2*tau2hat[1]*covar.sep(x, X3, d=1/b[1], g=g) + rho[2]^2*tau2hat[2]*covar.sep(x, X3, d=1/b[2], g=g) + tau2hat[3]*covar.sep(x, X3, d=1/b[3], g=g))
  
  V1 <- tau2hat[1]*covar.sep(X1, d=1/b[1], g=g)
  V12 <- rho[1]*tau2hat[1]*covar.sep(X1, X2, d=1/b[1], g=0)
  V13 <- rho[1]*rho[2]*tau2hat[1]*covar.sep(X1, X3, d=1/b[1], g=0)
  V2 <- rho[1]^2*tau2hat[1]*covar.sep(X2, d=1/b[1], g=g) + tau2hat[2]*covar.sep(X2, d=1/b[2], g=g)
  V23 <- rho[1]^2*rho[2]*tau2hat[1]*covar.sep(X2, X3, d=1/b[1], g=0) + rho[2]*tau2hat[2]*covar.sep(X2, X3, d=1/b[2], g=0)
  V3 <- rho[1]^2*rho[2]^2*tau2hat[1]*covar.sep(X3, d=1/b[1], g=g) + rho[2]^2*tau2hat[2]*covar.sep(X3, d=1/b[2], g=g) + tau2hat[3]*covar.sep(X3, d=1/b[3], g=g)
  
  V_3 <- rbind(cbind(V1, V12, V13), cbind(t(V12), V2, V23), cbind(t(V13), t(V23), V3))
  
  mx2 <- tx2 %*% solve(V_3) %*% c(Y1, Y2, Y3)
  
  koh.var2 <- pmax(0, diag(tau2hat[3]*covar.sep(matrix(x), d=1/b[3], g=g) + tau2hat[2]*rho[2]^2*covar.sep(matrix(x), d=1/b[2], g=g) + tau2hat[1]*rho[1]^2*rho[2]^2*covar.sep(matrix(x), d=1/b[3], g=g) - tx2 %*% solve(V_3+diag(g, nrow(V_3)))%*%t(tx2)))
  
  fit$X1 <- X1
  fit$Y1 <- Y1
  fit$X2 <- X2
  fit$Y2 <- Y2
  
  return(list(IMSPE=mean(koh.var2), fit=fit))
}


IMSPEKOH3select <- function(x, newx, fit, mc.sample=10000){
  ### This is updating KOH when the one data point is added ###
  # x; usually grid
  # newx; new point which will be added, should be n=1
  # fit; fitted KOH model
  X1 <- fit$X1
  X2 <- fit$X2
  X3 <- fit$X3
  Y1 <- fit$Y1
  Y2 <- fit$Y2
  Y3 <- fit$Y3
  
  b <- fit$b
  rho <- fit$rho
  tau2hat <- fit$tau2hat
  g <- fit$g
  
  ### Generate output
  y1.select <- fl(newx, l=1)
  
  X1 <- rbind(fit$X1, newx)
  Y1 <- rbind(fit$Y1, y1.select)  
  
  ### 2nd order KOH ###
  y2.select <- fl(newx, l=3)
  
  X2 <- rbind(fit$X2, newx)
  Y2 <- rbind(fit$Y2, y2.select)
  
  ### 3rd order KOH ###
  y3.select <- fl(newx, l=5)
  
  X3 <- rbind(fit$X3, newx)
  Y3 <- rbind(fit$Y3, y3.select)
  
  ### update sig2
  tx2 <- cbind(rho[1]*rho[2]*tau2hat[1]*covar.sep(x, X1, d=1/b[1], g=g), 
               rho[1]^2*rho[2]*tau2hat[1]*covar.sep(x, X2, d=1/b[1], g=g) + rho[2]*tau2hat[2]*covar.sep(x, X2, d=1/b[2], g=g),
               rho[1]^2*rho[2]^2*tau2hat[1]*covar.sep(x, X3, d=1/b[1], g=g) + rho[2]^2*tau2hat[2]*covar.sep(x, X3, d=1/b[2], g=g) + tau2hat[3]*covar.sep(x, X3, d=1/b[3], g=g))
  
  V1 <- tau2hat[1]*covar.sep(X1, d=1/b[1], g=g)
  V12 <- rho[1]*tau2hat[1]*covar.sep(X1, X2, d=1/b[1], g=0)
  V13 <- rho[1]*rho[2]*tau2hat[1]*covar.sep(X1, X3, d=1/b[1], g=0)
  V2 <- rho[1]^2*tau2hat[1]*covar.sep(X2, d=1/b[1], g=g) + tau2hat[2]*covar.sep(X2, d=1/b[2], g=g)
  V23 <- rho[1]^2*rho[2]*tau2hat[1]*covar.sep(X2, X3, d=1/b[1], g=0) + rho[2]*tau2hat[2]*covar.sep(X2, X3, d=1/b[2], g=0)
  V3 <- rho[1]^2*rho[2]^2*tau2hat[1]*covar.sep(X3, d=1/b[1], g=g) + rho[2]^2*tau2hat[2]*covar.sep(X3, d=1/b[2], g=g) + tau2hat[3]*covar.sep(X3, d=1/b[3], g=g)
  
  V_3 <- rbind(cbind(V1, V12, V13), cbind(t(V12), V2, V23), cbind(t(V13), t(V23), V3))
  
  mx2 <- tx2 %*% solve(V_3) %*% c(Y1, Y2, Y3)
  
  koh.var2 <- pmax(0, diag(tau2hat[3]*covar.sep(matrix(x), d=1/b[3], g=g) + tau2hat[2]*rho[2]^2*covar.sep(matrix(x), d=1/b[2], g=g) + tau2hat[1]*rho[1]^2*rho[2]^2*covar.sep(matrix(x), d=1/b[3], g=g) - tx2 %*% solve(V_3+diag(g, nrow(V_3)))%*%t(tx2)))
  
  fit$X1 <- X1
  fit$Y1 <- Y1  
  fit$X2 <- X2
  fit$Y2 <- Y2  
  fit$X3 <- X3
  fit$Y3 <- Y3
  
  return(list(IMSPE=mean(koh.var2), fit=fit))
}

