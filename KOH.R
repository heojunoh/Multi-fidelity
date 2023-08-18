library(plgp)
eps <- sqrt(.Machine$double.eps) 

KOHGP <- function(X, y, g=eps, center=TRUE){ # For 1st level's hyperparameter
  
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
  
  outg <- optim(rep(sort(distance(X))[sort(distance(X))!=0][0.1*length(sort(distance(X))[sort(distance(X))!=0])], ncol(X)), nlsep, method="L-BFGS-B", 
                lower=c(rep(0.001*sqrt(ncol(X)), ncol(X))), upper=c(rep(1000*sqrt(ncol(X)), ncol(X))), X=X, Y=y) 
  
  K <- covar.sep(X, d=outg$par, g=g)
  Ki <- solve(K+diag(g,n))
  tau2hat <- drop(t(y) %*% Ki %*% y / n)
  
  return(list(theta = outg$par, g=g, Ki=Ki, X = X, y = y, tau2hat=tau2hat))
}

KOH <- function(X, Y2, Y1, g=eps, center=TRUE){ # For after 1st level's hyperparameter
  
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
  
  outg <- optim(c(rep(sort(distance(X))[sort(distance(X))!=0][0.1*length(sort(distance(X))[sort(distance(X))!=0])], ncol(X)), 0.1), # For lengthscale and rho #
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

fit.KOH <- function(X1, X2, X3, Y1, Y2, Y3, g=eps){ # need to change function for another example
  
  ### KOH method ###
  Y2d3 <- fl(X3, l=3)
  Y1d2 <- fl(X2, l=1)
  
  ### estimating first order ###
  fit.KOHGP <- KOHGP(X1, Y1)
  b1 <- 1/fit.KOHGP$theta
  sig2_1 <- fit.KOHGP$tau2hat
  
  ### estimating second order ###
  KOHGP1 <- KOH(X2, Y2, Y1d2)
  rho1 <- KOHGP1$rho
  b2 <- 1/KOHGP1$theta
  sig2_2 <- KOHGP1$tau2hat
  
  ### estimating third order ###
  KOHGP2 <- KOH(X3, Y3, Y2d3)
  rho2 <- KOHGP2$rho
  b3 <- 1/KOHGP2$theta
  sig2_3 <- KOHGP2$tau2hat
  
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
  tx2 <- cbind(rho[1]*rho[2]*tau2hat[1]*covar.sep(x, X1, d=1/b[1], g=0), 
               rho[1]^2*rho[2]*tau2hat[1]*covar.sep(x, X2, d=1/b[1], g=0) + rho[2]*tau2hat[2]*covar.sep(x, X2, d=1/b[2], g=0),
               rho[1]^2*rho[2]^2*tau2hat[1]*covar.sep(x, X3, d=1/b[1], g=0) + rho[2]^2*tau2hat[2]*covar.sep(x, X3, d=1/b[2], g=0) + tau2hat[3]*covar.sep(x, X3, d=1/b[3], g=0))
  
  V1 <- tau2hat[1]*covar.sep(X1, d=1/b[1], g=g)
  V12 <- rho[1]*tau2hat[1]*covar.sep(X1, X2, d=1/b[1], g=0)
  V13 <- rho[1]*rho[2]*tau2hat[1]*covar.sep(X1, X3, d=1/b[1], g=0)
  V2 <- rho[1]^2*tau2hat[1]*covar.sep(X2, d=1/b[1], g=g) + tau2hat[2]*covar.sep(X2, d=1/b[2], g=g)
  V23 <- rho[1]^2*rho[2]*tau2hat[1]*covar.sep(X2, X3, d=1/b[1], g=0) + rho[2]*tau2hat[2]*covar.sep(X2, X3, d=1/b[2], g=0)
  V3 <- rho[1]^2*rho[2]^2*tau2hat[1]*covar.sep(X3, d=1/b[1], g=g) + rho[2]^2*tau2hat[2]*covar.sep(X3, d=1/b[2], g=g) + tau2hat[3]*covar.sep(X3, d=1/b[3], g=g)
  
  V_3 <- rbind(cbind(V1, V12, V13), cbind(t(V12), V2, V23), cbind(t(V13), t(V23), V3))
  
  mx2 <- tx2 %*% solve(V_3) %*% c(Y1, Y2, Y3)
  
  ### posterior variance ###
  koh.var2 <- pmax(0, diag(tau2hat[3]*covar.sep(matrix(x), d=1/b[3], g=0) + tau2hat[2]*rho[2]^2*covar.sep(matrix(x), d=1/b[2], g=0) + tau2hat[1]*rho[1]^2*rho[2]^2*covar.sep(matrix(x), d=1/b[3], g=0) - tx2 %*% solve(V_3+diag(g, nrow(V_3)))%*%t(tx2)))
  
  return(list(mu=mx2, sig2=koh.var2))
}


IMSPEKOH <- function(x, newx, fit, level){
  ### This is updating KOH when the one data point is added ###
  # x; usually grid
  # newx; new point which will be added, should be n=1
  # fit; fitted KOH model
  # level; level of fidelity
  X1 <- fit$X1
  X2 <- fit$X2
  
  b <- fit$b
  rho <- fit$rho
  tau2hat <- fit$tau2hat
  g <- fit$g
  
  
  if(level==1){ ### level 1 ###
    X1 <- rbind(fit$X1, newx)
    
    ### update sig2
    tx1 <- cbind(rho*tau2hat[1]*covar.sep(x, X1, d=1/b[1:2], g=g), 
                 rho^2*tau2hat[1]*covar.sep(x, X2, d=1/b[1:2], g=g) + tau2hat[2]*covar.sep(x, X2, d=1/b[3:4], g=g))
    V1 <- tau2hat[1]*covar.sep(X1, d=1/b[1:2], g=g)
    V12 <- rho*tau2hat[1]*covar.sep(X1, X2, d=1/b[1:2], g=0)
    V2 <- rho^2*tau2hat[1]*covar.sep(X2, d=1/b[1:2], g=g) + tau2hat[2]*covar.sep(X2, d=1/b[3:4], g=g)
    
    V_2 <- rbind(cbind(V1, V12), cbind(t(V12), V2))
    
    koh.var2 <- pmax(0, diag(tau2hat[2]*covar.sep(as.matrix(x), d=1/b[3:4], g=g) + tau2hat[1]*rho^2*covar.sep(as.matrix(x), d=1/b[1:2], g=g) - tx1 %*% solve(V_2)%*%t(tx1)))
  }else{ ### level 2 ###
    X1 <- rbind(fit$X1, newx)
    X2 <- rbind(fit$X2, newx)
    
    ### update sig2
    tx1 <- cbind(rho*tau2hat[1]*covar.sep(x, X1, d=1/b[1:2], g=g), 
                 rho^2*tau2hat[1]*covar.sep(x, X2, d=1/b[1:2], g=g) + tau2hat[2]*covar.sep(x, X2, d=1/b[3:4], g=g))
    V1 <- tau2hat[1]*covar.sep(X1, d=1/b[1:2], g=g)
    V12 <- rho*tau2hat[1]*covar.sep(X1, X2, d=1/b[1:2], g=0)
    V2 <- rho^2*tau2hat[1]*covar.sep(X2, d=1/b[1:2], g=g) + tau2hat[2]*covar.sep(X2, d=1/b[3:4], g=g)
    
    V_2 <- rbind(cbind(V1, V12), cbind(t(V12), V2))
    
    koh.var2 <- pmax(0, diag(tau2hat[2]*covar.sep(as.matrix(x), d=1/b[3:4], g=g) + tau2hat[1]*rho^2*covar.sep(as.matrix(x), d=1/b[1:2], g=g) - tx1 %*% solve(V_2)%*%t(tx1)))
  }
  
  return(IMSPE=mean(koh.var2))
}


IMSPEKOHselect <- function(x, newx, fit, level){
  ### This is updating KOH when the one data point is added ###
  # x; usually grid
  # newx; new point which will be added, should be n=1
  # fit; fitted KOH model
  X1 <- fit$X1
  X2 <- fit$X2
  Y1 <- fit$Y1
  Y2 <- fit$Y2
  
  newx <- matrix(newx, nrow=1)
  # b <- fit$b
  # rho <- fit$rho
  # tau2hat <- fit$tau2hat
  # g <- fit$g
  
  if(level==1){ ### level 1 ###
    ### Generate output
    d1 <- data.frame(newx*0.5+0.25, rep(0.05, 1)) # scale X to [-1,1]
    write.csv(d1, "/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/temp_to_matlab.txt", row.names=F)
    run_matlab_script("/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/SolveJetBlade.m", verbose = FALSE, desktop = FALSE, 
                      splash = FALSE, display = FALSE, wait = TRUE, single_thread = FALSE,
                      intern = TRUE)
    d2 <- read.table("/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/temp_to_r.txt", sep = ",")
    y1.select <- d2$V4
    
    X1 <- rbind(fit$X1, newx)
    Y1 <- c(fit$Y1, y1.select)
    
  }else{ ### level 2 ###
    ### Generate output
    d1 <- data.frame(newx*0.5+0.25, rep(0.05, 1)) # scale X to [-1,1]
    write.csv(d1, "/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/temp_to_matlab.txt", row.names=F)
    run_matlab_script("/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/SolveJetBlade.m", verbose = FALSE, desktop = FALSE, 
                      splash = FALSE, display = FALSE, wait = TRUE, single_thread = FALSE,
                      intern = TRUE)
    d2 <- read.table("/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/temp_to_r.txt", sep = ",")
    y1.select <- d2$V4
    
    d1 <- data.frame(newx*0.5+0.25, rep(0.025, 1)) # scale X to [-1,1]
    write.csv(d1, "/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/temp_to_matlab.txt", row.names=F)
    run_matlab_script("/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/SolveJetBlade.m", verbose = FALSE, desktop = FALSE, 
                      splash = FALSE, display = FALSE, wait = TRUE, single_thread = FALSE,
                      intern = TRUE)
    d2 <- read.table("/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/temp_to_r.txt", sep = ",")
    y2.select <- d2$V4
    
    X1 <- rbind(fit$X1, newx)
    X2 <- rbind(fit$X2, newx)
    Y1 <- c(fit$Y1, y1.select)
    Y2 <- c(fit$Y2, y2.select)
  }
  
  fit <- fit.KOH(X1, X2, Y1, Y2)
  
  return(fit=fit)
}


IMSPEKOH1 <- function(x, newx, fit, level){
  ### This is updating KOH when the one data point is added ###
  # x; usually grid
  # newx; new point which will be added, should be n=1
  # fit; fitted KOH model
  # level; level of fidelity
  X1 <- fit$X1
  X2 <- fit$X2
  
  b <- fit$b
  rho <- fit$rho
  tau2hat <- fit$tau2hat
  g <- fit$g
  
  
  if(level==1){ ### level 1 ###
    X1 <- rbind(fit$X1, newx)
    
    ### update sig2
    tx1 <- cbind(rho*tau2hat[1]*covar.sep(x, X1, d=1/b[1:8], g=g), 
                 rho^2*tau2hat[1]*covar.sep(x, X2, d=1/b[1:8], g=g) + tau2hat[2]*covar.sep(x, X2, d=1/b[9:16], g=g))
    
    V1 <- tau2hat[1]*covar.sep(X1, d=1/b[1:8], g=g)
    V12 <- rho*tau2hat[1]*covar.sep(X1, X2, d=1/b[1:8], g=0)
    V2 <- rho^2*tau2hat[1]*covar.sep(X2, d=1/b[1:8], g=g) + tau2hat[2]*covar.sep(X2, d=1/b[9:16], g=g)
    
    V_2 <- rbind(cbind(V1, V12), cbind(t(V12), V2))
    
    koh.var1 <- pmax(0, diag(tau2hat[2]*covar.sep(as.matrix(x), d=1/b[9:16], g=g) + tau2hat[1]*rho^2*covar.sep(as.matrix(x), d=1/b[1:8], g=g) - tx1 %*% solve(V_2)%*%t(tx1)))
    
  }else if(level==2){ ### level 2 ###
    X2 <- rbind(fit$X2, newx)
    
    ### update sig2
    tx1 <- cbind(rho*tau2hat[1]*covar.sep(x, X1, d=1/b[1:8], g=g), 
                 rho^2*tau2hat[1]*covar.sep(x, X2, d=1/b[1:8], g=g) + tau2hat[2]*covar.sep(x, X2, d=1/b[9:16], g=g))
    
    V1 <- tau2hat[1]*covar.sep(X1, d=1/b[1:8], g=g)
    V12 <- rho*tau2hat[1]*covar.sep(X1, X2, d=1/b[1:8], g=0)
    V2 <- rho^2*tau2hat[1]*covar.sep(X2, d=1/b[1:8], g=g) + tau2hat[2]*covar.sep(X2, d=1/b[9:16], g=g)
    
    V_2 <- rbind(cbind(V1, V12), cbind(t(V12), V2))
    
    koh.var1 <- pmax(0, diag(tau2hat[2]*covar.sep(as.matrix(x), d=1/b[9:16], g=g) + tau2hat[1]*rho^2*covar.sep(as.matrix(x), d=1/b[1:8], g=g) - tx1 %*% solve(V_2)%*%t(tx1)))
  }else{ 
    stop("level should be 1 or 2.")
  }
  
  return(IMSPE=mean(koh.var1))
}


IMSPEKOHselect1 <- function(x, newx, fit, level){
  ### This is updating KOH when the one data point is added ###
  # x; usually grid
  # newx; new point which will be added, should be n=1
  # fit; fitted KOH model
  X1 <- fit$X1
  X2 <- fit$X2
  Y1 <- fit$Y1
  Y2 <- fit$Y2
  
  # b <- fit$b
  # rho <- fit$rho
  # tau2hat <- fit$tau2hat
  # g <- fit$g
  
  if(level==1){ ### level 1 ###
    ### Generate output
    # y.select <- f1(newx)
    # y.select <- park91alc(newx)
    y.select <- apply(newx,1,outputlow.f)
    
    X1 <- rbind(fit$X1, newx)
    Y1 <- c(fit$Y1, y.select)
    
  }else if(level==2){ ### level 2 ###
    ### Generate output
    # y.select <- f2(newx)
    # y.select <- park91a(newx)
    y.select <- apply(newx,1,output.f)
    
    X2 <- rbind(fit$X2, newx)
    Y2 <- c(fit$Y2, y.select)

  }else{ 
    stop("level should be 1 or 2.")
  }
  
  fit <- fit.KOH(X1, X2, Y1, Y2)
  
  return(fit=fit)
}
