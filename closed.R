### Function to calculate the closed form ###
closed <- function(X1, y1, X2, y2, constant=FALSE){ # fit1; first layer, fit2; last layer's emulator such as f_M(X2, f1(x2))
  if(constant){
    fit1 <- GP(X1, y1, constant=TRUE)
    fit2 <- GP(cbind(X2, pred.GP(fit1, X2)$mu), y2, constant=TRUE) 
  }else{
    fit1 <- GP(X1, y1)
    fit2 <- GP(cbind(X2, pred.GP(fit1, X2)$mu), y2) 
  }
  return(list(fit1=fit1, fit2=fit2, constant=constant))
}

predclosed <- function(fit, x){
  constant <- fit$constant
  fit1 <- fit$fit1
  fit2 <- fit$fit2
  
  if(constant){
    d <- ncol(fit1$X)
    x <- matrix(x, ncol=d)
    x.mu <- pred.GP(fit1, x)$mu # mean of f1(u)
    sig2 <- pred.GP(fit1, x)$sig2 #*0
    
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
    sig2 <- sig2/attr(fit2$X,"scaled:scale")[d+1]^2
    
    # mean
    predy <- c(rep(0, nrow(x)))
    
    for(j in 1: nrow(x)){ # each test point
      predv <- c(rep(1,n))
      for(i in 1:n){ # each row of train set
        for(m in 1:d){ # dim of train set
          predv[i] <- predv[i] * exp(-(x[j,m]-X2[i,m])^2/theta[m]) # common components
        } # depends on kernel structure
        predv[i] <- predv[i] * 1/sqrt(1+2*sig2[j]/theta[d+1]) *
          exp(-(w1.x2[i]-x.mu[j])^2/(theta[d+1]+2*sig2[j]))
      }
      predy[j] <- mu2 + drop(predv%*%a)
    }
    
    # var
    predsig2 <- c(rep(0, nrow(x)))
    
    for(i in 1: nrow(x)){ # each test point
      mat <- matrix(1, n, n)
      for(k in 1:n){ # each row of train set
        for(l in 1:n){ # dim of train set
          for(m in 1:d){
            mat[k,l] <- mat[k,l] * exp(-((x[i,m]-X2[k,m])^2+(x[i,m]-X2[l,m])^2)/theta[m]) # common components
          }
          mat[k,l] <- mat[k,l] * #(a[k]*a[l] - tau2hat*Ci[k,l]) *
            1/sqrt(1+4*sig2[i]/theta[d+1]) * 
            exp(-((w1.x2[k]+w1.x2[l])/2-x.mu[i])^2/(theta[d+1]/2+2*sig2[i])) * 
            exp(-(w1.x2[k]-w1.x2[l])^2/(2*theta[d+1]))
        }}
      
      predsig2[i] <- pmax(0, tau2hat - (predy[i]-mu2)^2 + drop(t(a)%*%mat%*%a) - tau2hat*sum(diag(Ci%*%mat)))
    }
  }else{
    d <- ncol(fit1$X)
    x <- matrix(x, ncol=d)
    x.mu <- pred.GP(fit1, x)$mu # mean of f1(u)
    sig2 <- pred.GP(fit1, x)$sig2
    
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
    sig2 <- sig2/attr(fit2$X,"scaled:scale")[d+1]^2
    
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
          mat[k,l] <- mat[k,l] * #(a[k]*a[l] - tau2hat*Ci[k,l]) *
            1/sqrt(1+4*sig2[i]/theta[d+1]) * 
            exp(-((w1.x2[k]+w1.x2[l])/2-x.mu[i])^2/(theta[d+1]/2+2*sig2[i])) * 
            exp(-(w1.x2[k]-w1.x2[l])^2/(2*theta[d+1]))
        }}
      
      predsig2[i] <- pmax(0, tau2hat - predy[i]^2 + drop(t(a)%*%mat%*%a) - tau2hat*sum(diag(Ci%*%mat)))
    }
  }
  return(list(mu=predy, sig2=predsig2))
}


closed2 <- function(X1, y1, X2, y2, X3, y3, constant=FALSE){ # fit1; first layer, fit2; second layer f_M(X2, f1(X2)), fit3; last layer f_H(X3, f_M(X3, f1(X3)))
  if(constant){
    closed1 <- closed(X1, y1, X2, y2, constant=TRUE)
    
    fit1 <- closed1$fit1
    fit2 <- closed1$fit2
    fit3 <- GP(cbind(X3, pred.GP(fit2, cbind(X3, pred.GP(fit1, X3)$mu))$mu), y3, constant=TRUE)
  }else{
    closed1 <- closed(X1, y1, X2, y2)
    
    fit1 <- closed1$fit1
    fit2 <- closed1$fit2
    fit3 <- GP(cbind(X3, pred.GP(fit2, cbind(X3, pred.GP(fit1, X3)$mu))$mu), y3)
  }
  return(list(closed1=closed1, fit3=fit3, constant=constant))
}

predclosed2 <- function(fit, x){
  constant <- fit$constant
  closed1 <- fit$closed1
  fit1 <- closed1$fit1
  fit2 <- closed1$fit2
  fit3 <- fit$fit3
  
  if(constant){
    closed1 <- predclosed(closed1, x)
    x.mu <- closed1$mu
    sig2 <- closed1$sig2 
    
    d <- ncol(fit1$X)
    x <- matrix(x, ncol=d)
    
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
    sig2 <- sig2/attr(fit3$X,"scaled:scale")[d+1]^2
    
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
  }else{
    closed1 <- predclosed(closed1, x)
    x.mu <- closed1$mu
    sig2 <- closed1$sig2 
    
    d <- ncol(fit1$X)
    x <- matrix(x, ncol=d)
    
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
    sig2 <- sig2/attr(fit3$X,"scaled:scale")[d+1]^2
    
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
  }
  return(list(mu=predy, sig2=predsig2))
}

