### Function to calculate the closed form ###
closed <- function(x, fit1, fit2, constant=FALSE){ # fit1; first layer, fit2; last layer's emulator such as f_M(X2, f1(x2))
  if(constant){
    d <- ncol(fit1$X)
    x <- matrix(x, ncol=d)
    x.mu <- pred.GP(fit1, x)$mu # mean of f1(u)
    sig2 <- pred.GP(fit1, x)$sig2 # variance of f1(x)
    
    ### calculate the closed form ###
    X2 <- matrix(fit2$X[,-(d+1)], ncol=d)
    w1.x2 <- fit2$X[,d+1]
    y2 <- fit2$y
    n <- length(y2)
    theta <- fit2$theta
    tau2hat <- fit2$tau2hat
    mu2 <- fit2$mu.hat
    
    Ci <- fit2$Ki
    a <- Ci %*% (y2 - mu2)
    
    ### scale new inputs ###
    x <- t((t(x)-attr(fit2$X,"scaled:center")[1:d])/attr(fit2$X,"scaled:scale")[1:d])
    x.mu <- t((t(x.mu)-attr(fit2$X,"scaled:center")[d+1])/attr(fit2$X,"scaled:scale")[d+1])
    
    # mean
    predy <- c(rep(0, nrow(x)))
    
    for(j in 1: nrow(x)){
      predv <- c(rep(1,n))
      for(i in 1:n){
        for(m in 1:d){
          predv[i] <- predv[i] * exp(-(x[j,m]-X2[i,m])^2/theta[m])
        }
        predv[i] <- predv[i] * 1/sqrt(1+2*sig2[j]/theta[d+1]) *
          exp(-(w1.x2[i]-x.mu[j])^2/(theta[d+1]+2*sig2[j]))
      }
      predy[j] <- mu2 + drop(predv%*%a)
    }
    
    # var
    predsig2 <- c(rep(0, nrow(x)))
    
    for(i in 1: nrow(x)){
      mat <- matrix(1, n, n)
      for(k in 1:n){
        for(l in 1:n){
          for(m in 1:d){
            mat[k,l] <- mat[k,l] * exp(-((x[i,m]-X2[k,m])^2+(x[i,m]-X2[l,m])^2)/theta[m]) 
          }
          mat[k,l] <- mat[k,l] * (a[k]*a[l] - tau2hat*Ci[k,l]) *
            1/sqrt(1+4*sig2[i]/theta[d+1]) * 
            exp(-((w1.x2[k]+w1.x2[l])/2-x.mu[i])^2/(theta[d+1]/2+2*sig2[i])) * 
            exp(-(w1.x2[k]-w1.x2[l])^2/(2*theta[d+1]))
        }}
      
      predsig2[i] <- pmax(0, tau2hat - (predy[i]-mu2)^2 + sum(mat))
    }
    
    return(list(mu=predy, sig2=predsig2))
  }else{
    d <- ncol(fit1$X)
    x <- matrix(x, ncol=d)
    x.mu <- pred.GP(fit1, x)$mu # mean of f1(u)
    sig2 <- pred.GP(fit1, x)$sig2 # variance of f1(x)
    
    ### calculate the closed form ###
    X2 <- matrix(fit2$X[,-(d+1)], ncol=d)
    w1.x2 <- fit2$X[,d+1]
    y2 <- fit2$y
    n <- length(y2)
    theta <- fit2$theta
    tau2hat <- fit2$tau2hat
    
    Ci <- fit2$Ki
    a <- Ci %*% (y2 + attr(y2, "scaled:center"))
    
    ### scale new inputs ###
    x <- t((t(x)-attr(fit2$X,"scaled:center")[1:d])/attr(fit2$X,"scaled:scale")[1:d])
    x.mu <- t((t(x.mu)-attr(fit2$X,"scaled:center")[d+1])/attr(fit2$X,"scaled:scale")[d+1])
    
    # mean
    predy <- c(rep(0, nrow(x)))
    
    for(j in 1: nrow(x)){
      predv <- c(rep(1,n))
      for(i in 1:n){
        for(m in 1:d){
          predv[i] <- predv[i] * exp(-(x[j,m]-X2[i,m])^2/theta[m])
        }
        predv[i] <- predv[i] * 1/sqrt(1+2*sig2[j]/theta[d+1]) *
          exp(-(w1.x2[i]-x.mu[j])^2/(theta[d+1]+2*sig2[j]))
      }
      predy[j] <- drop(predv%*%a)
    }
    
    # var
    predsig2 <- c(rep(0, nrow(x)))
    
    for(i in 1: nrow(x)){
      mat <- matrix(1, n, n)
      for(k in 1:n){
        for(l in 1:n){
          for(m in 1:d){
            mat[k,l] <- mat[k,l] * exp(-((x[i,m]-X2[k,m])^2+(x[i,m]-X2[l,m])^2)/theta[m]) 
          }
          mat[k,l] <- mat[k,l] * (a[k]*a[l] - tau2hat*Ci[k,l]) *
            1/sqrt(1+4*sig2[i]/theta[d+1]) * 
            exp(-((w1.x2[k]+w1.x2[l])/2-x.mu[i])^2/(theta[d+1]/2+2*sig2[i])) * 
            exp(-(w1.x2[k]-w1.x2[l])^2/(2*theta[d+1]))
        }}
      
      predsig2[i] <- pmax(0, tau2hat - predy[i]^2 + sum(mat))
    }
    
    return(list(mu=predy, sig2=predsig2))
  }
}

closed2 <- function(x, fit1, fit2, fit3, constant=FALSE){ # fit1; first layer, fit2; second layer f_M(X2, f1(X2)), fit3; last layer f_H(X3, f_M(X3, f1(X3)))
  if(constant){
    d <- ncol(fit1$X)
    x <- matrix(x, ncol=d)
    
    ### second layer's output of test data 
    closed1 <- closed(x, fit1, fit2, constant=TRUE)
    x.mu <- closed1$mu
    sig2 <- closed1$sig2
    
    ### calculate the closed form ###
    X3 <- matrix(fit3$X[,-(d+1)], ncol=d)
    w2.x3 <- fit3$X[,d+1]
    y3 <- fit3$y
    n <- length(y3)
    theta <- fit3$theta
    tau2hat <- fit3$tau2hat
    mu3 <- fit3$mu.hat
    
    Ci <- fit3$Ki
    a <- Ci %*% (y3 - mu3)
    
    ### scale new inputs ###
    x <- t((t(x)-attr(fit3$X,"scaled:center")[1:d])/attr(fit3$X,"scaled:scale")[1:d])
    x.mu <- t((t(x.mu)-attr(fit3$X,"scaled:center")[d+1])/attr(fit3$X,"scaled:scale")[d+1])
    
    # mean
    predy <- c(rep(0, nrow(x)))
    
    for(j in 1: nrow(x)){
      predv <- c(rep(1,n))
      for(i in 1:n){
        for(m in 1:d){
          predv[i] <- predv[i] * exp(-(x[j,m]-X3[i,m])^2/theta[m])
        }
        predv[i] <- predv[i] * 1/sqrt(1+2*sig2[j]/theta[d+1]) *
          exp(-(w2.x3[i]-x.mu[j])^2/(theta[d+1]+2*sig2[j]))
      }
      predy[j] <- mu3 + drop(predv%*%a)
    }
    
    # var
    predsig2 <- c(rep(0, nrow(x)))
    
    for(i in 1: nrow(x)){
      mat <- matrix(1, n, n)
      for(k in 1:n){
        for(l in 1:n){
          for(m in 1:d){
            mat[k,l] <- mat[k,l] * exp(-((x[i,m]-X3[k,m])^2+(x[i,m]-X3[l,m])^2)/theta[m]) 
          }
          mat[k,l] <- mat[k,l] * (a[k]*a[l] - tau2hat*Ci[k,l]) *
            1/sqrt(1+4*sig2[i]/theta[d+1]) * 
            exp(-((w2.x3[k]+w2.x3[l])/2-x.mu[i])^2/(theta[d+1]/2+2*sig2[i])) * 
            exp(-(w2.x3[k]-w2.x3[l])^2/(2*theta[d+1]))
        }}
      
      predsig2[i] <- pmax(0, tau2hat - (predy[i]-mu3)^2 + sum(mat))
    }
    
    return(list(mu=predy, sig2=predsig2))
  }else{
    d <- ncol(fit1$X)
    x <- matrix(x, ncol=d)
    
    ### second layer's output of test data 
    closed1 <- closed(x, fit1, fit2)
    x.mu <- closed1$mu
    sig2 <- closed1$sig2
    
    ### calculate the closed form ###
    X3 <- matrix(fit3$X[,-(d+1)], ncol=d)
    w2.x3 <- fit3$X[,d+1]
    y3 <- fit3$y
    n <- length(y3)
    theta <- fit3$theta
    tau2hat <- fit3$tau2hat
    
    Ci <- fit3$Ki
    a <- Ci %*% (y3 + attr(y3, "scaled:center"))
    
    ### scale new inputs ###
    x <- t((t(x)-attr(fit3$X,"scaled:center")[1:d])/attr(fit3$X,"scaled:scale")[1:d])
    x.mu <- t((t(x.mu)-attr(fit3$X,"scaled:center")[d+1])/attr(fit3$X,"scaled:scale")[d+1])
    
    # mean
    predy <- c(rep(0, nrow(x)))
    
    for(j in 1: nrow(x)){
      predv <- c(rep(1,n))
      for(i in 1:n){
        for(m in 1:d){
          predv[i] <- predv[i] * exp(-(x[j,m]-X3[i,m])^2/theta[m])
        }
        predv[i] <- predv[i] * 1/sqrt(1+2*sig2[j]/theta[d+1]) *
          exp(-(w2.x3[i]-x.mu[j])^2/(theta[d+1]+2*sig2[j]))
      }
      predy[j] <- drop(predv%*%a)
    }
    
    # var
    predsig2 <- c(rep(0, nrow(x)))
    
    for(i in 1: nrow(x)){
      mat <- matrix(1, n, n)
      for(k in 1:n){
        for(l in 1:n){
          for(m in 1:d){
            mat[k,l] <- mat[k,l] * exp(-((x[i,m]-X3[k,m])^2+(x[i,m]-X3[l,m])^2)/theta[m]) 
          }
          mat[k,l] <- mat[k,l] * (a[k]*a[l] - tau2hat*Ci[k,l]) *
            1/sqrt(1+4*sig2[i]/theta[d+1]) * 
            exp(-((w2.x3[k]+w2.x3[l])/2-x.mu[i])^2/(theta[d+1]/2+2*sig2[i])) * 
            exp(-(w2.x3[k]-w2.x3[l])^2/(2*theta[d+1]))
        }}
      
      predsig2[i] <- pmax(0, tau2hat - predy[i]^2 + sum(mat))
    }
    
    return(list(mu=predy, sig2=predsig2))
  }
}
