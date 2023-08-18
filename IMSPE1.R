### Function to calculate IMSPE of n+1 when the low fidelity is added ###
library(lhs)
library(laGP)
source("GP.R")
source("closed.R")

### For closed1 ###

IMSPEfast1 <- function(x, newx, fit, mc.sample=10000, level){ 
  ### This is updating Ki using Ki_{n+1}'s closed form ###
  # x; usually grid
  # newx; new point which will be added, should be n=1
  # fit1; first layer's emulator
  # fit2; second layer's emulator f_M(X2, f1(x2))
  
  fit1 <- fit$fit1
  fit2 <- fit$fit2
  constant <- fit$constant
  newx <- matrix(newx, nrow=1)
  
  ### Choose level 1 ###
  if(level == 1){
    x1.sample <- rnorm(mc.sample, mean=pred.GP(fit1, newx)$mu, sd=sqrt(pred.GP(fit1, newx)$sig2))
    mu.cand1 <- mean(x1.sample) # f1(newx)
    
    x.center1 <- attr(fit1$X, "scaled:center")
    x.scale1 <- attr(fit1$X, "scaled:scale")
    y.center1 <- attr(fit1$y, "scaled:center")
    
    newx <- matrix((newx-attr(fit1$X,"scaled:center"))/attr(fit1$X,"scaled:scale"), nrow=1)
    
    ### update Ki
    v.next <- drop(covar.sep(X1=newx, d=fit1$theta, g=0) - 
                     t(covar.sep(X1=fit1$X, X2=newx, d=fit1$theta, g=0)) %*%
                     fit1$Ki %*%
                     covar.sep(X1=fit1$X, X2=newx, d=fit1$theta, g=0))
    g.next <- - drop(solve(v.next)) * fit1$Ki %*% covar.sep(X1=fit1$X, X2=newx, d=fit1$theta, g=0)
    fit1$Ki <- rbind(cbind(fit1$Ki+g.next%*%t(g.next)*v.next, g.next),
                     cbind(t(g.next), solve(v.next)))
    
    fit1$X <- rbind(fit1$X, newx)
    attr(fit1$X, "scaled:center") <- x.center1
    attr(fit1$X, "scaled:scale") <- x.scale1
    
    if(constant){
      fit1$y <- c(fit1$y, mu.cand1)
    }else{
      fit1$y <- c(fit1$y, mu.cand1-attr(fit1$y,"scaled:center"))
      attr(fit1$y, "scaled:center") <- y.center1
    }
    
    fit1$tau2hat <- drop(t(fit1$y) %*% fit1$Ki %*% fit1$y / length(fit1$y))
    
    fit$fit1 <- fit1
    fit$fit2 <- fit2
    
    if(constant){
      return(list(IMSPE=mean(predclosed(fit, x)$sig2), fit1new=fit1))
    }else{
      return(list(IMSPE=mean(predclosed(fit, x)$sig2), fit1new=fit1))
    }
  }
  

  ### Choose level 2 ###
  if(level == 2){
    x1.sample <- rnorm(mc.sample, mean=pred.GP(fit1, newx)$mu, sd=sqrt(pred.GP(fit1, newx)$sig2))
    mu.cand1 <- mean(x1.sample) # f1(newx)
    
    x2.sample <- rnorm(mc.sample, 
                       mean=pred.GP(fit2, cbind(newx, mu.cand1))$mu, 
                       sd=sqrt(pred.GP(fit2, cbind(newx, mu.cand1))$sig2))
    mu.cand2 <- mean(x2.sample) # f2(newx)
    
    x.center1 <- attr(fit1$X, "scaled:center")
    x.scale1 <- attr(fit1$X, "scaled:scale")
    y.center1 <- attr(fit1$y, "scaled:center")
    
    x.center2 <- attr(fit2$X, "scaled:center")
    x.scale2 <- attr(fit2$X, "scaled:scale")
    y.center2 <- attr(fit2$y, "scaled:center")
    
    newx1 <- matrix((newx-attr(fit1$X,"scaled:center"))/attr(fit1$X,"scaled:scale"), nrow=1) 
    newx2 <- t((t(cbind(newx, mean(x1.sample)))-attr(fit2$X,"scaled:center"))/attr(fit2$X,"scaled:scale"))
    
    ### update Ki1
    v.next1 <- drop(covar.sep(X1=newx1, d=fit1$theta, g=0) - 
                      t(covar.sep(X1=fit1$X, X2=newx1, d=fit1$theta, g=0)) %*%
                      fit1$Ki %*%
                      covar.sep(X1=fit1$X, X2=newx1, d=fit1$theta, g=0))
    g.next1 <- - drop(solve(v.next1)) * fit1$Ki %*% covar.sep(X1=fit1$X, X2=newx1, d=fit1$theta, g=0)
    fit1$Ki <- rbind(cbind(fit1$Ki+g.next1%*%t(g.next1)*v.next1, g.next1),
                     cbind(t(g.next1), solve(v.next1)))
    
    fit1$X <- rbind(fit1$X, newx1)
    attr(fit1$X, "scaled:center") <- x.center1
    attr(fit1$X, "scaled:scale") <- x.scale1
    
    if(constant){
      fit1$y <- c(fit1$y, mu.cand1)
    }else{
      fit1$y <- c(fit1$y, mu.cand1-attr(fit1$y,"scaled:center"))
      attr(fit1$y, "scaled:center") <- y.center1
    }
    
    fit1$tau2hat <- drop(t(fit1$y) %*% fit1$Ki %*% fit1$y / length(fit1$y))
    
    ### update Ki2
    v.next2 <- drop(covar.sep(X1=newx2, d=fit2$theta, g=0) - 
                      t(covar.sep(X1=fit2$X, X2=newx2, d=fit2$theta, g=0)) %*%
                      fit2$Ki %*%
                      covar.sep(X1=fit2$X, X2=newx2, d=fit2$theta, g=0))
    g.next2 <- - drop(solve(v.next2)) * fit2$Ki %*% covar.sep(X1=fit2$X, X2=newx2, d=fit2$theta, g=0)
    fit2$Ki <- rbind(cbind(fit2$Ki+g.next2%*%t(g.next2)*v.next2, g.next2),
                     cbind(t(g.next2), solve(v.next2)))
    
    fit2$X <- rbind(fit2$X, newx2)
    attr(fit2$X, "scaled:center") <- x.center2
    attr(fit2$X, "scaled:scale") <- x.scale2
    
    if(constant){
      fit2$y <- c(fit2$y, mu.cand2)
    }else{
      fit2$y <- c(fit2$y, mu.cand2-attr(fit2$y,"scaled:center"))
      attr(fit2$y, "scaled:center") <- y.center2
    }
    
    fit2$tau2hat <- drop(t(fit2$y) %*% fit2$Ki %*% fit2$y / length(fit2$y))
    
    fit$fit1 <- fit1
    fit$fit2 <- fit2
    
    if(constant){
      return(list(IMSPE=mean(predclosed(fit, x)$sig2), fit1new=fit1, fit2new=fit2))
    }else{
      return(list(IMSPE=mean(predclosed(fit, x)$sig2), fit1new=fit1, fit2new=fit2))
    }
  }
}


IMSPEselect1 <- function(x, newx, closedfit, mc.sample=10000, level){ 
  ### Updating Ki using true output (not MC sample) ###
  # x; usually grid
  # newx; new point which will be added, should be n=1
  # fit1; first layer's emulator
  # fit2; second layer's emulator f_M(X2, f1(x2))
  
  fit1 <- closedfit$fit1
  fit2 <- closedfit$fit2
  constant <- closedfit$constant
  
  newx <- matrix(newx, nrow=1)
  
  X1 <- t(t(fit1$X)*attr(fit1$X,"scaled:scale")+attr(fit1$X,"scaled:center"))
  X2 <- matrix(t(t(fit2$X)*attr(fit2$X,"scaled:scale")+attr(fit2$X,"scaled:center"))[,-ncol(fit2$X)], ncol=ncol(fit2$X)-1)
  
  if(constant){
    y1 <- fit1$y
    y2 <- fit2$y
  }else{
    y1 <- fit1$y+attr(fit1$y,"scaled:center")
    y2 <- fit2$y+attr(fit2$y,"scaled:center")
  }
  
  
  ### Choose level 1 ###
  if(level == 1){
    y1.select <- f1(newx) ### need to change for different function
    # y1.select <- park91alc(newx)
    
    X1 <- rbind(X1, newx)
    y1 <- c(y1, y1.select)
  }
  
  ### Choose level 2 ###
  if(level == 2){
    y1.select <- f1(newx) ### need to change for different function
    y2.select <- f2(newx)
    # y1.select <- park91alc(newx)
    # y2.select <- park91a(newx)
    
    X1 <- rbind(X1, newx)
    y1 <- c(y1, y1.select)
    X2 <- rbind(X2, newx)
    y2 <- c(y2, y2.select)
  }
  
  fit <- closed(X1, y1, X2, y2, constant)
  
  return(list(IMSPE=mean(predclosed(fit, x)$sig2), fit=fit))
}


### For closed2 ###

IMSPEfast2 <- function(x, newx, fit, mc.sample=10000, level){ 
  ### This is updating Ki using Ki_{n+1}'s closed form ###
  # x; usually grid
  # newx; new point which will be added, should be n=1 or nrow=1
  # fit1; first layer's emulator
  # fit2; second layer's emulator f_M(X2, f1(x2))
  
  fit1 <- fit$fit1
  fit2 <- fit$fit2
  constant <- fit$constant
  newx <- matrix(newx, nrow=1)
  
  ### Choose level 1 ###
  if(level == 1){
    x1.sample <- rnorm(mc.sample, mean=pred.GP(fit1, newx)$mu, sd=sqrt(pred.GP(fit1, newx)$sig2))
    mu.cand1 <- mean(x1.sample) # f1(newx)
    
    x.center1 <- attr(fit1$X, "scaled:center")
    x.scale1 <- attr(fit1$X, "scaled:scale")
    y.center1 <- attr(fit1$y, "scaled:center")
    
    newx <- matrix((newx-attr(fit1$X,"scaled:center"))/attr(fit1$X,"scaled:scale"), nrow=1)
    
    ### update Ki
    v.next <- drop(covar.sep(X1=newx, d=fit1$theta, g=0) - 
                     t(covar.sep(X1=fit1$X, X2=newx, d=fit1$theta, g=0)) %*%
                     fit1$Ki %*%
                     covar.sep(X1=fit1$X, X2=newx, d=fit1$theta, g=0))
    g.next <- - drop(solve(v.next)) * fit1$Ki %*% covar.sep(X1=fit1$X, X2=newx, d=fit1$theta, g=0)
    fit1$Ki <- rbind(cbind(fit1$Ki+g.next%*%t(g.next)*v.next, g.next),
                     cbind(t(g.next), solve(v.next)))
    
    fit1$X <- rbind(fit1$X, newx)
    attr(fit1$X, "scaled:center") <- x.center1
    attr(fit1$X, "scaled:scale") <- x.scale1
    
    if(constant){
      fit1$y <- c(fit1$y, mu.cand1)
    }else{
      fit1$y <- c(fit1$y, mu.cand1-attr(fit1$y,"scaled:center"))
      attr(fit1$y, "scaled:center") <- y.center1
    }
    
    fit1$tau2hat <- drop(t(fit1$y-fit1$mu.hat) %*% fit1$Ki %*% (fit1$y-fit1$mu.hat) / length(fit1$y))
    
    fit$fit1 <- fit1
    fit$fit2 <- fit2
    
    if(constant){
      return(list(IMSPE=mean(predclosed(fit, x)$sig2), fit1new=fit1))
    }else{
      return(list(IMSPE=mean(predclosed(fit, x)$sig2), fit1new=fit1))
    }
  }
  
  
  ### Choose level 2 ###
  if(level == 2){
    x1.sample <- rnorm(mc.sample, mean=pred.GP(fit1, newx)$mu, sd=sqrt(pred.GP(fit1, newx)$sig2))
    mu.cand1 <- mean(x1.sample) # f1(newx)
    
    x2.sample <- rnorm(mc.sample, 
                       mean=pred.GP(fit2, cbind(newx, mu.cand1))$mu, 
                       sd=sqrt(pred.GP(fit2, cbind(newx, mu.cand1))$sig2))
    mu.cand2 <- mean(x2.sample) # f2(newx)
    
    x.center1 <- attr(fit1$X, "scaled:center")
    x.scale1 <- attr(fit1$X, "scaled:scale")
    y.center1 <- attr(fit1$y, "scaled:center")
    
    x.center2 <- attr(fit2$X, "scaled:center")
    x.scale2 <- attr(fit2$X, "scaled:scale")
    y.center2 <- attr(fit2$y, "scaled:center")
    
    newx1 <- matrix((newx-attr(fit1$X,"scaled:center"))/attr(fit1$X,"scaled:scale"), nrow=1) 
    newx2 <- t((t(cbind(newx, mean(x1.sample)))-attr(fit2$X,"scaled:center"))/attr(fit2$X,"scaled:scale"))
    
    ### update Ki1
    v.next1 <- drop(covar.sep(X1=newx1, d=fit1$theta, g=0) - 
                      t(covar.sep(X1=fit1$X, X2=newx1, d=fit1$theta, g=0)) %*%
                      fit1$Ki %*%
                      covar.sep(X1=fit1$X, X2=newx1, d=fit1$theta, g=0))
    g.next1 <- - drop(solve(v.next1)) * fit1$Ki %*% covar.sep(X1=fit1$X, X2=newx1, d=fit1$theta, g=0)
    fit1$Ki <- rbind(cbind(fit1$Ki+g.next1%*%t(g.next1)*v.next1, g.next1),
                     cbind(t(g.next1), solve(v.next1)))
    
    fit1$X <- rbind(fit1$X, newx1)
    attr(fit1$X, "scaled:center") <- x.center1
    attr(fit1$X, "scaled:scale") <- x.scale1
    
    if(constant){
      fit1$y <- c(fit1$y, mu.cand1)
    }else{
      fit1$y <- c(fit1$y, mu.cand1-attr(fit1$y,"scaled:center"))
      attr(fit1$y, "scaled:center") <- y.center1
    }
    
    fit1$tau2hat <- drop(t(fit1$y-fit1$mu.hat) %*% fit1$Ki %*% (fit1$y-fit1$mu.hat) / length(fit1$y))
    
    ### update Ki2
    v.next2 <- drop(covar.sep(X1=newx2, d=fit2$theta, g=0) - 
                      t(covar.sep(X1=fit2$X, X2=newx2, d=fit2$theta, g=0)) %*%
                      fit2$Ki %*%
                      covar.sep(X1=fit2$X, X2=newx2, d=fit2$theta, g=0))
    g.next2 <- - drop(solve(v.next2)) * fit2$Ki %*% covar.sep(X1=fit2$X, X2=newx2, d=fit2$theta, g=0)
    fit2$Ki <- rbind(cbind(fit2$Ki+g.next2%*%t(g.next2)*v.next2, g.next2),
                     cbind(t(g.next2), solve(v.next2)))
    
    fit2$X <- rbind(fit2$X, newx2)
    attr(fit2$X, "scaled:center") <- x.center2
    attr(fit2$X, "scaled:scale") <- x.scale2
    
    if(constant){
      fit2$y <- c(fit2$y, mu.cand2)
    }else{
      fit2$y <- c(fit2$y, mu.cand2-attr(fit2$y,"scaled:center"))
      attr(fit2$y, "scaled:center") <- y.center2
    }
    
    fit2$tau2hat <- drop(t(fit2$y-fit2$mu.hat) %*% fit2$Ki %*% (fit2$y-fit2$mu.hat) / length(fit2$y))
    
    fit$fit1 <- fit1
    fit$fit2 <- fit2
    
    if(constant){
      return(list(IMSPE=mean(predclosed(fit, x)$sig2), fit1new=fit1, fit2new=fit2))
    }else{
      return(list(IMSPE=mean(predclosed(fit, x)$sig2), fit1new=fit1, fit2new=fit2))
    }
  }
}


IMSPEselect2 <- function(x, newx, closedfit, mc.sample=10000, level){ 
  ### Updating Ki using true output (not MC sample) ###
  # x; usually grid
  # newx; new point which will be added, should be n=1
  # fit1; first layer's emulator
  # fit2; second layer's emulator f_M(X2, f1(x2))
  
  fit1 <- closedfit$fit1
  fit2 <- closedfit$fit2
  constant <- closedfit$constant
  
  newx <- matrix(newx, nrow=1)
  
  X1 <- t(t(fit1$X)*attr(fit1$X,"scaled:scale")+attr(fit1$X,"scaled:center"))
  X2 <- matrix(t(t(fit2$X)*attr(fit2$X,"scaled:scale")+attr(fit2$X,"scaled:center"))[,-ncol(fit2$X)], ncol=ncol(fit2$X)-1)
  
  if(constant){
    y1 <- fit1$y
    y2 <- fit2$y
  }else{
    y1 <- fit1$y+attr(fit1$y,"scaled:center")
    y2 <- fit2$y+attr(fit2$y,"scaled:center")
  }
  
  
  ### Choose level 1 ###
  if(level == 1){
    # y1.select <- f1(newx) ### need to change for different function
    d1 <- data.frame(newx*0.5+0.25, rep(0.05, 1)) # scale X to [-1,1]
    write.csv(d1, "/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/temp_to_matlab.txt", row.names=F)
    run_matlab_script("/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/SolveJetBlade.m", verbose = FALSE, desktop = FALSE, 
                      splash = FALSE, display = FALSE, wait = TRUE, single_thread = FALSE,
                      intern = TRUE)
    d2 <- read.table("/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/temp_to_r.txt", sep = ",")
    y1.select <- d2$V4
    
    X1 <- rbind(X1, newx)
    y1 <- c(y1, y1.select)
  }
  
  ### Choose level 2 ###
  if(level == 2){
    # y1.select <- f1(newx) ### need to change for different function
    # y2.select <- f2(newx) 
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
    
    X1 <- rbind(X1, newx)
    y1 <- c(y1, y1.select)
    X2 <- rbind(X2, newx)
    y2 <- c(y2, y2.select)
  }
  
  fit <- closed(X1, y1, X2, y2, constant)
  
  return(list(IMSPE=mean(predclosed(fit, x)$sig2), fit=fit))
}


