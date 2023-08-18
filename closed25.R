### Function to calculate the closed form ###
closed25 <- function(X1, y1, X2, y2, constant=FALSE){ # fit1; first layer, fit2; last layer's emulator such as f_M(X2, f1(x2))
  if(constant){
    fit1 <- matGP(X1, y1, nu=2.5, constant=TRUE)
    fit2 <- matGP(cbind(X2, pred.matGP(fit1, X2)$mu), y2, nu=2.5, constant=TRUE) 
  }else{
    fit1 <- matGP(X1, y1, nu=2.5)
    fit2 <- matGP(cbind(X2, pred.matGP(fit1, X2)$mu), y2, nu=2.5) 
  }
  return(list(fit1=fit1, fit2=fit2, constant=constant))
}


predclosed25 <- function(fit, x){
  constant <- fit$constant
  fit1 <- fit$fit1
  fit2 <- fit$fit2
  
  if(constant){
    d <- ncol(fit1$X)
    x <- matrix(x, ncol=d)
    x.mu <- pred.matGP(fit1, x)$mu # mean of f1(u)
    sig2 <- pred.matGP(fit1, x)$sig2 #*0
    
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
          predv[i] <- predv[i] * matern.kernel(sqrt(distance(t(t(x[j,m])/theta[m]), t(t(X2[i,m])/theta[m]))), nu=2.5) # common but depends on kernel
        } # depends on kernel structure
        
        mua <- x.mu[j] - sqrt(5)*sig2[j]/theta[d+1]
        mub <- x.mu[j] + sqrt(5)*sig2[j]/theta[d+1]
        
        lambda11 <- c(1, mua, mua^2+sig2[j])
        lambda12 <- c(0, 1, mua+w1.x2[i])
        lambda21 <- c(1, -mub, mub^2+sig2[j])
        lambda22 <- c(0, 1, -mub-w1.x2[i])
        
        e1 <- c(1-sqrt(5)*w1.x2[i]/theta[d+1]+5*w1.x2[i]^2/(3*theta[d+1]^2),
                sqrt(5)/theta[d+1]-10*w1.x2[i]/(3*theta[d+1]^2),
                5/(3*theta[d+1]^2))
        e2 <- c(1+sqrt(5)*w1.x2[i]/theta[d+1]+5*w1.x2[i]^2/(3*theta[d+1]^2),
                sqrt(5)/theta[d+1]+10*w1.x2[i]/(3*theta[d+1]^2),
                5/(3*theta[d+1]^2))
        
        predv[i] <- predv[i] * (exp((5*sig2[j] + 2*sqrt(5)*theta[d+1]*(w1.x2[i] - x.mu[j]))/(2*theta[d+1]^2)) *
                                  (e1 %*% lambda11 * pnorm((mua - w1.x2[i])/sqrt(sig2[j])) +
                                     e1 %*% lambda12 * sqrt(sig2[j])/sqrt(2*pi) * exp(-(w1.x2[i] - mua)^2/(2*sig2[j]))) +
                                  exp((5*sig2[j] - 2*sqrt(5)*theta[d+1]*(w1.x2[i] - x.mu[j]))/(2*theta[d+1]^2)) *
                                  (e2 %*% lambda21 * pnorm((-mub + w1.x2[i])/sqrt(sig2[j])) +
                                     e2 %*% lambda22 * sqrt(sig2[j])/sqrt(2*pi) * exp(-(w1.x2[i] - mub)^2/(2*sig2[j]))))
      }
      predy[j] <- mu2 + drop(predv%*%a)
    }
    
    # var
    predsig2 <- c(rep(0, nrow(x)))
    
    for(i in 1: nrow(x)){
      mat <- matrix(1, n, n)
      for(k in 1:n){
        for(l in 1:n){
          for(m in 1:d){ # common but depends on kernel
            mat[k,l] <- mat[k,l] * matern.kernel(sqrt(distance(x[i,m]/theta[m], X2[k,m]/theta[m])), nu=2.5) *
              matern.kernel(sqrt(distance(x[i,m]/theta[m], X2[l,m]/theta[m])), nu=2.5) 
          } # expected depends on kernel structure
          mat[k,l] <- mat[k,l] * zetafun(w1=w1.x2[k], w2=w1.x2[l], m=x.mu[i], s=sig2[i], nu=2.5, theta=theta[d+1])
        }}
      
      predsig2[i] <- pmax(0, tau2hat - (predy[i]-mu2)^2 + drop(t(a)%*%mat%*%a) - tau2hat*sum(diag(Ci%*%mat)))
    }
  }else{
    d <- ncol(fit1$X)
    x <- matrix(x, ncol=d)
    x.mu <- pred.matGP(fit1, x)$mu # mean of f1(u)
    sig2 <- pred.matGP(fit1, x)$sig2
    
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
    
    for(j in 1: nrow(x)){ # each test point
      predv <- c(rep(1,n))
      for(i in 1:n){ # each row of train set
        for(m in 1:d){ # dim of train set
          predv[i] <- predv[i] * matern.kernel(sqrt(distance(t(t(x[j,m])/theta[m]), t(t(X2[i,m])/theta[m]))), nu=2.5) # common but depends on kernel
        } # depends on kernel structure
        
        mua <- x.mu[j] - sqrt(5)*sig2[j]/theta[d+1]
        mub <- x.mu[j] + sqrt(5)*sig2[j]/theta[d+1]
        
        lambda11 <- c(1, mua, mua^2+sig2[j])
        lambda12 <- c(0, 1, mua+w1.x2[i])
        lambda21 <- c(1, -mub, mub^2+sig2[j])
        lambda22 <- c(0, 1, -mub-w1.x2[i])
        
        e1 <- c(1-sqrt(5)*w1.x2[i]/theta[d+1]+5*w1.x2[i]^2/(3*theta[d+1]^2),
                sqrt(5)/theta[d+1]-10*w1.x2[i]/(3*theta[d+1]^2),
                5/(3*theta[d+1]^2))
        e2 <- c(1+sqrt(5)*w1.x2[i]/theta[d+1]+5*w1.x2[i]^2/(3*theta[d+1]^2),
                sqrt(5)/theta[d+1]+10*w1.x2[i]/(3*theta[d+1]^2),
                5/(3*theta[d+1]^2))
        
        predv[i] <- predv[i] * (exp((5*sig2[j] + 2*sqrt(5)*theta[d+1]*(w1.x2[i] - x.mu[j]))/(2*theta[d+1]^2)) *
                                  (e1 %*% lambda11 * pnorm((mua - w1.x2[i])/sqrt(sig2[j])) +
                                     e1 %*% lambda12 * sqrt(sig2[j])/sqrt(2*pi) * exp(-(w1.x2[i] - mua)^2/(2*sig2[j]))) +
                                  exp((5*sig2[j] - 2*sqrt(5)*theta[d+1]*(w1.x2[i] - x.mu[j]))/(2*theta[d+1]^2)) *
                                  (e2 %*% lambda21 * pnorm((-mub + w1.x2[i])/sqrt(sig2[j])) +
                                     e2 %*% lambda22 * sqrt(sig2[j])/sqrt(2*pi) * exp(-(w1.x2[i] - mub)^2/(2*sig2[j]))))
      }
      predy[j] <- drop(predv%*%a)
    }
    
    # var
    predsig2 <- c(rep(0, nrow(x)))
    
    for(i in 1: nrow(x)){
      mat <- matrix(1, n, n)
      for(k in 1:n){
        for(l in 1:n){
          for(m in 1:d){ # common but depends on kernel
            mat[k,l] <- mat[k,l] * matern.kernel(sqrt(distance(t(t(x[i,m])/theta[m]), t(t(X2[k,m])/theta[m]))), nu=2.5) *
              matern.kernel(sqrt(distance(t(t(x[i,m])/theta[m]), t(t(X2[l,m])/theta[m]))), nu=2.5) 
          } # expected depends on kernel structure
          mat[k,l] <- mat[k,l] * zetafun(w1=w1.x2[k], w2=w1.x2[l], m=x.mu[i], s=sig2[i], nu=2.5, theta=theta[d+1])
        }}
      
      predsig2[i] <- pmax(0, tau2hat - predy[i]^2 + drop(t(a)%*%mat%*%a) - tau2hat*sum(diag(Ci%*%mat)))
    }
  }
  return(list(mu=predy, sig2=predsig2))
}


closed225 <- function(X1, y1, X2, y2, X3, y3, constant=FALSE){ # fit1; first layer, fit2; second layer f_M(X2, f1(X2)), fit3; last layer f_H(X3, f_M(X3, f1(X3)))
  if(constant){
    closed1 <- closed25(X1, y1, X2, y2, constant=TRUE)
    
    fit1 <- closed1$fit1
    fit2 <- closed1$fit2
    fit3 <- matGP(cbind(X3, pred.matGP(fit2, cbind(X3, pred.matGP(fit1, X3)$mu))$mu), y3, nu=2.5, constant=TRUE)
  }else{
    closed1 <- closed25(X1, y1, X2, y2)
    
    fit1 <- closed1$fit1
    fit2 <- closed1$fit2
    fit3 <- matGP(cbind(X3, pred.matGP(fit2, cbind(X3, pred.matGP(fit1, X3)$mu))$mu), y3, nu=2.5)
  }
  return(list(closed1=closed1, fit3=fit3, constant=constant))
}


predclosed225 <- function(fit, x){
  constant <- fit$constant
  closed1 <- fit$closed1
  fit1 <- closed1$fit1
  fit2 <- closed1$fit2
  fit3 <- fit$fit3
  
  if(constant){
    d <- ncol(fit1$X)
    x <- matrix(x, ncol=d)
    closed1 <- predclosed25(closed1, x)
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
    sig2 <- sig2/attr(fit3$X,"scaled:scale")[d+1]^2
    
    # mean
    predy <- c(rep(0, nrow(x)))
    
    for(j in 1: nrow(x)){
      predv <- c(rep(1,n))
      for(i in 1:n){
        for(m in 1:d){
          predv[i] <- predv[i] * matern.kernel(sqrt(distance(t(t(x[j,m])/theta[m]), t(t(X3[i,m])/theta[m]))), nu=2.5) # common but depends on kernel
        } # depends on kernel structure

        mua <- x.mu[j] - sqrt(5)*sig2[j]/theta[d+1]
        mub <- x.mu[j] + sqrt(5)*sig2[j]/theta[d+1]
        
        lambda11 <- c(1, mua, mua^2+sig2[j])
        lambda12 <- c(0, 1, mua+w2.x3[i])
        lambda21 <- c(1, -mub, mub^2+sig2[j])
        lambda22 <- c(0, 1, -mub-w2.x3[i])
        
        e1 <- c(1-sqrt(5)*w2.x3[i]/theta[d+1]+5*w2.x3[i]^2/(3*theta[d+1]^2),
                sqrt(5)/theta[d+1]-10*w2.x3[i]/(3*theta[d+1]^2),
                5/(3*theta[d+1]^2))
        e2 <- c(1+sqrt(5)*w2.x3[i]/theta[d+1]+5*w2.x3[i]^2/(3*theta[d+1]^2),
                sqrt(5)/theta[d+1]+10*w2.x3[i]/(3*theta[d+1]^2),
                5/(3*theta[d+1]^2))
        
        predv[i] <- predv[i] * (exp((5*sig2[j] + 2*sqrt(5)*theta[d+1]*(w2.x3[i] - x.mu[j]))/(2*theta[d+1]^2)) *
                                  (e1 %*% lambda11 * pnorm((mua - w2.x3[i])/sqrt(sig2[j])) +
                                     e1 %*% lambda12 * sqrt(sig2[j])/sqrt(2*pi) * exp(-(w2.x3[i] - mua)^2/(2*sig2[j]))) +
                                  exp((5*sig2[j] - 2*sqrt(5)*theta[d+1]*(w2.x3[i] - x.mu[j]))/(2*theta[d+1]^2)) *
                                  (e2 %*% lambda21 * pnorm((-mub + w2.x3[i])/sqrt(sig2[j])) +
                                     e2 %*% lambda22 * sqrt(sig2[j])/sqrt(2*pi) * exp(-(w2.x3[i] - mub)^2/(2*sig2[j]))))
      }
      predy[j] <- mu3 + drop(predv%*%a)
    }
    
    # var
    predsig2 <- c(rep(0, nrow(x)))
    
    for(i in 1: nrow(x)){
      mat <- matrix(1, n, n)
      for(k in 1:n){
        for(l in 1:n){
          for(m in 1:d){ # common but depends on kernel
            mat[k,l] <- mat[k,l] * matern.kernel(sqrt(distance(t(t(x[i,m])/theta[m]), t(t(X3[k,m])/theta[m]))), nu=2.5) *
              matern.kernel(sqrt(distance(t(t(x[i,m])/theta[m]), t(t(X3[l,m])/theta[m]))), nu=2.5) 
          } # expected depends on kernel structure
          mat[k,l] <- mat[k,l] * zetafun(w1=w2.x3[k], w2=w2.x3[l], m=x.mu[i], s=sig2[i], nu=2.5, theta=theta[d+1])
        }}
      
      predsig2[i] <- pmax(0, tau2hat - (predy[i]-mu3)^2 + drop(t(a)%*%mat%*%a) - tau2hat*sum(diag(Ci%*%mat)))
    }
  }else{
    d <- ncol(fit1$X)
    x <- matrix(x, ncol=d)
    closed1 <- predclosed25(closed1, x)
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
    sig2 <- sig2/attr(fit3$X,"scaled:scale")[d+1]^2
    
    # mean
    predy <- c(rep(0, nrow(x)))
    
    for(j in 1: nrow(x)){
      predv <- c(rep(1,n))
      for(i in 1:n){
        for(m in 1:d){
          predv[i] <- predv[i] * matern.kernel(sqrt(distance(t(t(x[j,m])/theta[m]), t(t(X3[i,m])/theta[m]))), nu=2.5) # common but depends on kernel
        } # depends on kernel structure
        
        mua <- x.mu[j] - sqrt(5)*sig2[j]/theta[d+1]
        mub <- x.mu[j] + sqrt(5)*sig2[j]/theta[d+1]
        
        lambda11 <- c(1, mua, mua^2+sig2[j])
        lambda12 <- c(0, 1, mua+w2.x3[i])
        lambda21 <- c(1, -mub, mub^2+sig2[j])
        lambda22 <- c(0, 1, -mub-w2.x3[i])
        
        e1 <- c(1-sqrt(5)*w2.x3[i]/theta[d+1]+5*w2.x3[i]^2/(3*theta[d+1]^2),
                sqrt(5)/theta[d+1]-10*w2.x3[i]/(3*theta[d+1]^2),
                5/(3*theta[d+1]^2))
        e2 <- c(1+sqrt(5)*w2.x3[i]/theta[d+1]+5*w2.x3[i]^2/(3*theta[d+1]^2),
                sqrt(5)/theta[d+1]+10*w2.x3[i]/(3*theta[d+1]^2),
                5/(3*theta[d+1]^2))
        
        predv[i] <- predv[i] * (exp((5*sig2[j] + 2*sqrt(5)*theta[d+1]*(w2.x3[i] - x.mu[j]))/(2*theta[d+1]^2)) *
                                  (e1 %*% lambda11 * pnorm((mua - w2.x3[i])/sqrt(sig2[j])) +
                                     e1 %*% lambda12 * sqrt(sig2[j])/sqrt(2*pi) * exp(-(w2.x3[i] - mua)^2/(2*sig2[j]))) +
                                  exp((5*sig2[j] - 2*sqrt(5)*theta[d+1]*(w2.x3[i] - x.mu[j]))/(2*theta[d+1]^2)) *
                                  (e2 %*% lambda21 * pnorm((-mub + w2.x3[i])/sqrt(sig2[j])) +
                                     e2 %*% lambda22 * sqrt(sig2[j])/sqrt(2*pi) * exp(-(w2.x3[i] - mub)^2/(2*sig2[j]))))
      }
      predy[j] <- drop(predv%*%a)
    }
    
    # var
    predsig2 <- c(rep(0, nrow(x)))
    
    for(i in 1: nrow(x)){
      mat <- matrix(1, n, n)
      for(k in 1:n){
        for(l in 1:n){
          for(m in 1:d){ # common but depends on kernel
            mat[k,l] <- mat[k,l] * matern.kernel(sqrt(distance(t(t(x[i,m])/theta[m]), t(t(X3[k,m])/theta[m]))), nu=2.5) *
              matern.kernel(sqrt(distance(t(t(x[i,m])/theta[m]), t(t(X3[l,m])/theta[m]))), nu=2.5) 
          } # expected depends on kernel structure
          mat[k,l] <- mat[k,l] * zetafun(w1=w2.x3[k], w2=w2.x3[l], m=x.mu[i], s=sig2[i], nu=2.5, theta=theta[d+1])
        }}
      
      predsig2[i] <- pmax(0, tau2hat - predy[i]^2 + drop(t(a)%*%mat%*%a) - tau2hat*sum(diag(Ci%*%mat)))
    }
  }
  return(list(mu=predy, sig2=predsig2))
}


zetafun <- function(w1, w2, m, s, nu, theta){
  if(nu == 1.5){
    if(w1 < w2){
      muc <- m - 2*sqrt(3)*s/theta
      mud <- m + 2*sqrt(3)*s/theta
      
      lambda31 <- c(1, muc, muc^2+s)
      lambda32 <- c(0, 1, muc+w2)
      lambda41 <- c(1, m, m^2+s)
      lambda42 <- c(0, 1, m+w1)
      lambda43 <- c(0, 1, m+w2)
      lambda51 <- c(1, -mud, mud^2+s)
      lambda52 <- c(0, 1, -mud-w1)
      
      e3 <- c(1+(3*w1*w2 - sqrt(3)*theta*(w1+w2))/theta^2,
              (2*sqrt(3)*theta - 3*(w1+w2))/theta^2,
              3/theta^2)
      e4 <- c(1+(-3*w1*w2 + sqrt(3)*theta*(w2-w1))/theta^2,
              3*(w1+w2)/theta^2,
              -3/theta^2)
      e5 <- c(1+(3*w1*w2 + sqrt(3)*theta*(w1+w2))/theta^2,
              (2*sqrt(3)*theta + 3*(w1+w2))/theta^2,
              3/theta^2)
      
      
      return(exp((6*s + sqrt(3)*theta*(w1 + w2 - 2*m))/theta^2) *
               (e3 %*% lambda31 * pnorm((muc - w2)/sqrt(s)) +
                  e3 %*% lambda32 * sqrt(s)/sqrt(2*pi) * exp(-(w2 - muc)^2/(2*s))) +
               exp(-(sqrt(3)*(w2-w1)/theta)) *
               (e4 %*% lambda41 * (pnorm((w2 - m)/sqrt(s)) - pnorm((w1 - m)/sqrt(s))) +
                  e4 %*% lambda42 * sqrt(s)/sqrt(2*pi) * exp(-(w1 - m)^2/(2*s)) -
                  e4 %*% lambda43 * sqrt(s)/sqrt(2*pi) * exp(-(w2 - m)^2/(2*s)) ) +
               exp((6*s - sqrt(3)*theta*(w1 + w2 - 2*m))/theta^2) *
               (e5 %*% lambda51 * pnorm((w1 - mud)/sqrt(s)) +
                  e5 %*% lambda52 * sqrt(s)/sqrt(2*pi) * exp(-(w1 - mud)^2/(2*s))))
    }else{ # b <= a
      muc <- m - 2*sqrt(3)*s/theta
      mud <- m + 2*sqrt(3)*s/theta
      
      lambda31 <- c(1, muc, muc^2+s)
      lambda32 <- c(0, 1, muc+w1)
      lambda41 <- c(1, m, m^2+s)
      lambda42 <- c(0, 1, m+w2)
      lambda43 <- c(0, 1, m+w1)
      lambda51 <- c(1, -mud, mud^2+s)
      lambda52 <- c(0, 1, -mud-w2)
      
      e3 <- c(1+(3*w1*w2 - sqrt(3)*theta*(w1+w2))/theta^2,
              (2*sqrt(3)*theta - 3*(w1+w2))/theta^2,
              3/theta^2)
      e4 <- c(1+(-3*w1*w2 + sqrt(3)*theta*(w1-w2))/theta^2,
              3*(w1+w2)/theta^2,
              -3/theta^2)
      e5 <- c(1+(3*w1*w2 + sqrt(3)*theta*(w1+w2))/theta^2,
              (2*sqrt(3)*theta + 3*(w1+w2))/theta^2,
              3/theta^2)
      
      
      return(exp((6*s + sqrt(3)*theta*(w1 + w2 - 2*m))/theta^2) *
               (e3 %*% lambda31 * pnorm((muc - w1)/sqrt(s)) +
                  e3 %*% lambda32 * sqrt(s)/sqrt(2*pi) * exp(-(w1 - muc)^2/(2*s))) +
               exp(-(sqrt(3)*(w1-w2)/theta)) *
               (e4 %*% lambda41 * (pnorm((w1 - m)/sqrt(s)) - pnorm((w2 - m)/sqrt(s))) +
                  e4 %*% lambda42 * sqrt(s)/sqrt(2*pi) * exp(-(w2 - m)^2/(2*s)) -
                  e4 %*% lambda43 * sqrt(s)/sqrt(2*pi) * exp(-(w1 - m)^2/(2*s)) ) +
               exp((6*s - sqrt(3)*theta*(w1 + w2 - 2*m))/theta^2) *
               (e5 %*% lambda51 * pnorm((w2 - mud)/sqrt(s)) +
                  e5 %*% lambda52 * sqrt(s)/sqrt(2*pi) * exp(-(w2 - mud)^2/(2*s))))
    }
  }else if(nu == 2.5){
    if(w1 < w2){
      muc <- m - 2*sqrt(5)*s/theta
      mud <- m + 2*sqrt(5)*s/theta
      
      lambda31 <- c(1, muc, muc^2+s, muc^3+3*muc*s, muc^4+6*muc^2*s+3*s^2)
      lambda32 <- c(0, 1, muc+w2, muc^2+2*s+w2^2+muc*w2, muc^3+w2^3+muc^2*w2+muc*w2^2+3*s*w2+5*muc*s)
      lambda41 <- c(1, m, m^2+s, m^3+3*m*s, m^4+6*m^2*s+3*s^2)
      lambda42 <- c(0, 1, m+w1, m^2+2*s+w1^2+m*w1, m^3+w1^3+w1*m^2+m*w1^2+3*s*w1+5*s*m)
      lambda43 <- c(0, 1, m+w2, m^2+2*s+w2^2+m*w2, m^3+w2^3+w2*m^2+m*w2^2+3*s*w2+5*s*m)
      lambda51 <- c(1, -mud, mud^2+s, -mud^3-3*mud*s, mud^4+6*mud^2*s+3*s^2)
      lambda52 <- c(0, 1, -mud-w1, mud^2+2*s+w1^2+mud*w1, -mud^3-w1^3-mud^2*w1-mud*w1^2-3*s*w1-5*mud*s)
      
      e3 <- c(1+(25*w1^2*w2^2 - 3*sqrt(5)*(3*theta^3+5*theta*w1*w2)*(w1+w2) + 15*theta^2*(w1^2+w2^2+3*w1*w2))/(9*theta^4),
              (18*sqrt(5)*theta^3 + 15*sqrt(5)*theta*(w1^2+w2^2) - 75*theta^2*(w1+w2) - 50*w1*w2*(w1+w2) + 60*sqrt(5)*theta*w1*w2)/(9*theta^4),
              5*(5*w1^2+5*w2^2 + 15*theta^2 - 9*sqrt(5)*theta*(w1+w2) + 20*w1*w2)/(9*theta^4),
              10*(3*sqrt(5)*theta - 5*(w1+w2))/(9*theta^4),
              25/(9*theta^4))
      e4 <- c(1+(25*w1^2*w2^2 + 3*sqrt(5)*(3*theta^3-5*theta*w1*w2)*(w2-w1) + 15*theta^2*(w1^2+w2^2-3*w1*w2))/(9*theta^4),
              5*(3*sqrt(5)*theta*(w2^2-w1^2) + 3*theta^2*(w1+w2) - 10*w1*w2*(w1+w2))/(9*theta^4),
              5*(5*w1^2+5*w2^2 - 3*theta^2 - 3*sqrt(5)*theta*(w2-w1) + 20*w1*w2)/(9*theta^4),
              -50*(w1+w2)/(9*theta^4),
              25/(9*theta^4))
      e5 <- c(1+(25*w1^2*w2^2 + 3*sqrt(5)*(3*theta^3+5*theta*w1*w2)*(w1+w2) + 15*theta^2*(w1^2+w2^2+3*w1*w2))/(9*theta^4),
              (18*sqrt(5)*theta^3 + 15*sqrt(5)*theta*(w1^2+w2^2) + 75*theta^2*(w1+w2) + 50*w1*w2*(w1+w2) + 60*sqrt(5)*theta*w1*w2)/(9*theta^4),
              5*(5*w1^2+5*w2^2 + 15*theta^2 + 9*sqrt(5)*theta*(w1+w2) + 20*w1*w2)/(9*theta^4),
              10*(3*sqrt(5)*theta + 5*(w1+w2))/(9*theta^4),
              25/(9*theta^4))
      
      
      return(exp((10*s + sqrt(5)*theta*(w1 + w2 - 2*m))/theta^2) *
               (e3 %*% lambda31 * pnorm((muc - w2)/sqrt(s)) +
                  e3 %*% lambda32 * sqrt(s)/sqrt(2*pi) * exp(-(w2 - muc)^2/(2*s))) +
               exp(-(sqrt(5)*(w2-w1)/theta)) *
               (e4 %*% lambda41 * (pnorm((w2 - m)/sqrt(s)) - pnorm((w1 - m)/sqrt(s))) +
                  e4 %*% lambda42 * sqrt(s)/sqrt(2*pi) * exp(-(w1 - m)^2/(2*s)) -
                  e4 %*% lambda43 * sqrt(s)/sqrt(2*pi) * exp(-(w2 - m)^2/(2*s)) ) +
               exp((10*s - sqrt(5)*theta*(w1 + w2 - 2*m))/theta^2) *
               (e5 %*% lambda51 * pnorm((w1 - mud)/sqrt(s)) +
                  e5 %*% lambda52 * sqrt(s)/sqrt(2*pi) * exp(-(w1 - mud)^2/(2*s))))
    }else{
      muc <- m - 2*sqrt(5)*s/theta
      mud <- m + 2*sqrt(5)*s/theta
      
      lambda31 <- c(1, muc, muc^2+s, muc^3+3*muc*s, muc^4+6*muc^2*s+3*s^2)
      lambda32 <- c(0, 1, muc+w1, muc^2+2*s+w1^2+muc*w1, muc^3+w1^3+muc^2*w1+muc*w1^2+3*s*w1+5*muc*s)
      lambda41 <- c(1, m, m^2+s, m^3+3*m*s, m^4+6*m^2*s+3*s^2)
      lambda42 <- c(0, 1, m+w2, m^2+2*s+w2^2+m*w2, m^3+w2^3+w2*m^2+m*w2^2+3*s*w2+5*s*m)
      lambda43 <- c(0, 1, m+w1, m^2+2*s+w1^2+m*w1, m^3+w1^3+w1*m^2+m*w1^2+3*s*w1+5*s*m)
      lambda51 <- c(1, -mud, mud^2+s, -mud^3-3*mud*s, mud^4+6*mud^2*s+3*s^2)
      lambda52 <- c(0, 1, -mud-w2, mud^2+2*s+w2^2+mud*w2, -mud^3-w2^3-mud^2*w2-mud*w2^2-3*s*w2-5*mud*s)
      
      e3 <- c(1+(25*w1^2*w2^2 - 3*sqrt(5)*(3*theta^3+5*theta*w1*w2)*(w1+w2) + 15*theta^2*(w1^2+w2^2+3*w1*w2))/(9*theta^4),
              (18*sqrt(5)*theta^3 + 15*sqrt(5)*theta*(w1^2+w2^2) - 75*theta^2*(w1+w2) - 50*w1*w2*(w1+w2) + 60*sqrt(5)*theta*w1*w2)/(9*theta^4),
              5*(5*w1^2+5*w2^2 + 15*theta^2 - 9*sqrt(5)*theta*(w1+w2) + 20*w1*w2)/(9*theta^4),
              10*(3*sqrt(5)*theta - 5*(w1+w2))/(9*theta^4),
              25/(9*theta^4))
      e4 <- c(1+(25*w1^2*w2^2 + 3*sqrt(5)*(3*theta^3-5*theta*w1*w2)*(w1-w2) + 15*theta^2*(w1^2+w2^2-3*w1*w2))/(9*theta^4),
              5*(3*sqrt(5)*theta*(w1^2-w2^2) + 3*theta^2*(w1+w2) - 10*w1*w2*(w1+w2))/(9*theta^4),
              5*(5*w1^2+5*w2^2 - 3*theta^2 - 3*sqrt(5)*theta*(w1-w2) + 20*w1*w2)/(9*theta^4),
              -50*(w1+w2)/(9*theta^4),
              25/(9*theta^4))
      e5 <- c(1+(25*w1^2*w2^2 + 3*sqrt(5)*(3*theta^3+5*theta*w1*w2)*(w1+w2) + 15*theta^2*(w1^2+w2^2+3*w1*w2))/(9*theta^4),
              (18*sqrt(5)*theta^3 + 15*sqrt(5)*theta*(w1^2+w2^2) + 75*theta^2*(w1+w2) + 50*w1*w2*(w1+w2) + 60*sqrt(5)*theta*w1*w2)/(9*theta^4),
              5*(5*w1^2+5*w2^2 + 15*theta^2 + 9*sqrt(5)*theta*(w1+w2) + 20*w1*w2)/(9*theta^4),
              10*(3*sqrt(5)*theta + 5*(w1+w2))/(9*theta^4),
              25/(9*theta^4))
      
      
      return(exp((10*s + sqrt(5)*theta*(w1 + w2 - 2*m))/theta^2) *
               (e3 %*% lambda31 * pnorm((muc - w1)/sqrt(s)) +
                  e3 %*% lambda32 * sqrt(s)/sqrt(2*pi) * exp(-(w1 - muc)^2/(2*s))) +
               exp(-(sqrt(5)*(w1-w2)/theta)) *
               (e4 %*% lambda41 * (pnorm((w1 - m)/sqrt(s)) - pnorm((w2 - m)/sqrt(s))) +
                  e4 %*% lambda42 * sqrt(s)/sqrt(2*pi) * exp(-(w2 - m)^2/(2*s)) -
                  e4 %*% lambda43 * sqrt(s)/sqrt(2*pi) * exp(-(w1 - m)^2/(2*s)) ) +
               exp((10*s - sqrt(5)*theta*(w1 + w2 - 2*m))/theta^2) *
               (e5 %*% lambda51 * pnorm((w2 - mud)/sqrt(s)) +
                  e5 %*% lambda52 * sqrt(s)/sqrt(2*pi) * exp(-(w2 - mud)^2/(2*s))))
    }
  }
}


# zetafun25 <- function(w1, w2, m, s){ # i < j -> w1.x2[k] < w1.x2[l]
#   # w1=w1.x2[k], w2=w1.x2[l], m=x.mu[i], s=sig2[i]
#   if(w1 < w2){
#     muc <- m - 2*sqrt(5)*s/theta[d+1]
#     mud <- m + 2*sqrt(5)*s/theta[d+1]
#     
#     lambda31 <- c(1, muc, muc^2+s, muc^3+3*muc*s, muc^4+6*muc^2*s+3*s^2)
#     lambda32 <- c(0, 1, muc+w2, muc^2+2*s+w2^2+muc*w2, muc^3+w2^3+muc^2*w2+muc*w2^2+3*s*w2+5*muc*s)
#     lambda41 <- c(1, m, m^2+s, m^3+3*m*s, m^4+6*m^2*s+3*s^2)
#     lambda42 <- c(0, 1, m+w1, m^2+2*s+w1^2+m*w1, m^3+w1^3+w1*m^2+m*w1^2+3*s*w1+5*s*m)
#     lambda43 <- c(0, 1, m+w2, m^2+2*s+w2^2+m*w2, m^3+w2^3+w2*m^2+m*w2^2+3*s*w2+5*s*m)
#     lambda51 <- c(1, -mud, mud^2+s, -mud^3-3*mud*s, mud^4+6*mud^2*s+3*s^2)
#     lambda52 <- c(0, 1, -mud-w1, mud^2+2*s+w1^2+mud*w1, -mud^3-w1^3-mud^2*w1-mud*w1^2-3*s*w1-5*mud*s)
#     
#     e3 <- c(1+(25*w1^2*w2^2 - 3*sqrt(5)*(3*theta[d+1]^3+5*theta[d+1]*w1*w2)*(w1+w2) + 15*theta[d+1]^2*(w1^2+w2^2+3*w1*w2))/(9*theta[d+1]^4),
#             (18*sqrt(5)*theta[d+1]^3 + 15*sqrt(5)*theta[d+1]*(w1^2+w2^2) - 75*theta[d+1]^2*(w1+w2) - 50*w1*w2*(w1+w2) + 60*sqrt(5)*theta[d+1]*w1*w2)/(9*theta[d+1]^4),
#             5*(5*w1^2+5*w2^2 + 15*theta[d+1]^2 - 9*sqrt(5)*theta[d+1]*(w1+w2) + 20*w1*w2)/(9*theta[d+1]^4),
#             10*(3*sqrt(5)*theta[d+1] - 5*(w1+w2))/(9*theta[d+1]^4),
#             25/(9*theta[d+1]^4))
#     e4 <- c(1+(25*w1^2*w2^2 + 3*sqrt(5)*(3*theta[d+1]^3-5*theta[d+1]*w1*w2)*(w2-w1) + 15*theta[d+1]^2*(w1^2+w2^2-3*w1*w2))/(9*theta[d+1]^4),
#             5*(3*sqrt(5)*theta[d+1]*(w2^2-w1^2) + 3*theta[d+1]^2*(w1+w2) - 10*w1*w2*(w1+w2))/(9*theta[d+1]^4),
#             5*(5*w1^2+5*w2^2 - 3*theta[d+1]^2 - 3*sqrt(5)*theta[d+1]*(w2-w1) + 20*w1*w2)/(9*theta[d+1]^4),
#             -50*(w1+w2)/(9*theta[d+1]^4),
#             25/(9*theta[d+1]^4))
#     e5 <- c(1+(25*w1^2*w2^2 + 3*sqrt(5)*(3*theta[d+1]^3+5*theta[d+1]*w1*w2)*(w1+w2) + 15*theta[d+1]^2*(w1^2+w2^2+3*w1*w2))/(9*theta[d+1]^4),
#             (18*sqrt(5)*theta[d+1]^3 + 15*sqrt(5)*theta[d+1]*(w1^2+w2^2) + 75*theta[d+1]^2*(w1+w2) + 50*w1*w2*(w1+w2) + 60*sqrt(5)*theta[d+1]*w1*w2)/(9*theta[d+1]^4),
#             5*(5*w1^2+5*w2^2 + 15*theta[d+1]^2 + 9*sqrt(5)*theta[d+1]*(w1+w2) + 20*w1*w2)/(9*theta[d+1]^4),
#             10*(3*sqrt(5)*theta[d+1] + 5*(w1+w2))/(9*theta[d+1]^4),
#             25/(9*theta[d+1]^4))
#     
#     
#     return(exp((10*s + sqrt(5)*theta[d+1]*(w1 + w2 - 2*m))/theta[d+1]^2) *
#              (e3 %*% lambda31 * pnorm((muc - w2)/sqrt(s)) +
#                 e3 %*% lambda32 * sqrt(s)/sqrt(2*pi) * exp(-(w2 - muc)^2/(2*s))) +
#              exp(-(sqrt(5)*(w2-w1)/theta[d+1])) *
#              (e4 %*% lambda41 * (pnorm((w2 - m)/sqrt(s)) - pnorm((w1 - m)/sqrt(s))) +
#                 e4 %*% lambda42 * sqrt(s)/sqrt(2*pi) * exp(-(w1 - m)^2/(2*s)) -
#                 e4 %*% lambda43 * sqrt(s)/sqrt(2*pi) * exp(-(w2 - m)^2/(2*s)) ) +
#              exp((10*s - sqrt(5)*theta[d+1]*(w1 + w2 - 2*m))/theta[d+1]^2) *
#              (e5 %*% lambda51 * pnorm((w1 - mud)/sqrt(s)) +
#                 e5 %*% lambda52 * sqrt(s)/sqrt(2*pi) * exp(-(w1 - mud)^2/(2*s))))
#   }else{
#     muc <- m - 2*sqrt(5)*s/theta[d+1]
#     mud <- m + 2*sqrt(5)*s/theta[d+1]
#     
#     lambda31 <- c(1, muc, muc^2+s, muc^3+3*muc*s, muc^4+6*muc^2*s+3*s^2)
#     lambda32 <- c(0, 1, muc+w1, muc^2+2*s+w1^2+muc*w1, muc^3+w1^3+muc^2*w1+muc*w1^2+3*s*w1+5*muc*s)
#     lambda41 <- c(1, m, m^2+s, m^3+3*m*s, m^4+6*m^2*s+3*s^2)
#     lambda42 <- c(0, 1, m+w2, m^2+2*s+w2^2+m*w2, m^3+w2^3+w2*m^2+m*w2^2+3*s*w2+5*s*m)
#     lambda43 <- c(0, 1, m+w1, m^2+2*s+w1^2+m*w1, m^3+w1^3+w1*m^2+m*w1^2+3*s*w1+5*s*m)
#     lambda51 <- c(1, -mud, mud^2+s, -mud^3-3*mud*s, mud^4+6*mud^2*s+3*s^2)
#     lambda52 <- c(0, 1, -mud-w2, mud^2+2*s+w2^2+mud*w2, -mud^3-w2^3-mud^2*w2-mud*w2^2-3*s*w2-5*mud*s)
#     
#     e3 <- c(1+(25*w1^2*w2^2 - 3*sqrt(5)*(3*theta[d+1]^3+5*theta[d+1]*w1*w2)*(w1+w2) + 15*theta[d+1]^2*(w1^2+w2^2+3*w1*w2))/(9*theta[d+1]^4),
#             (18*sqrt(5)*theta[d+1]^3 + 15*sqrt(5)*theta[d+1]*(w1^2+w2^2) - 75*theta[d+1]^2*(w1+w2) - 50*w1*w2*(w1+w2) + 60*sqrt(5)*theta[d+1]*w1*w2)/(9*theta[d+1]^4),
#             5*(5*w1^2+5*w2^2 + 15*theta[d+1]^2 - 9*sqrt(5)*theta[d+1]*(w1+w2) + 20*w1*w2)/(9*theta[d+1]^4),
#             10*(3*sqrt(5)*theta[d+1] - 5*(w1+w2))/(9*theta[d+1]^4),
#             25/(9*theta[d+1]^4))
#     e4 <- c(1+(25*w1^2*w2^2 + 3*sqrt(5)*(3*theta[d+1]^3-5*theta[d+1]*w1*w2)*(w1-w2) + 15*theta[d+1]^2*(w1^2+w2^2-3*w1*w2))/(9*theta[d+1]^4),
#             5*(3*sqrt(5)*theta[d+1]*(w1^2-w2^2) + 3*theta[d+1]^2*(w1+w2) - 10*w1*w2*(w1+w2))/(9*theta[d+1]^4),
#             5*(5*w1^2+5*w2^2 - 3*theta[d+1]^2 - 3*sqrt(5)*theta[d+1]*(w1-w2) + 20*w1*w2)/(9*theta[d+1]^4),
#             -50*(w1+w2)/(9*theta[d+1]^4),
#             25/(9*theta[d+1]^4))
#     e5 <- c(1+(25*w1^2*w2^2 + 3*sqrt(5)*(3*theta[d+1]^3+5*theta[d+1]*w1*w2)*(w1+w2) + 15*theta[d+1]^2*(w1^2+w2^2+3*w1*w2))/(9*theta[d+1]^4),
#             (18*sqrt(5)*theta[d+1]^3 + 15*sqrt(5)*theta[d+1]*(w1^2+w2^2) + 75*theta[d+1]^2*(w1+w2) + 50*w1*w2*(w1+w2) + 60*sqrt(5)*theta[d+1]*w1*w2)/(9*theta[d+1]^4),
#             5*(5*w1^2+5*w2^2 + 15*theta[d+1]^2 + 9*sqrt(5)*theta[d+1]*(w1+w2) + 20*w1*w2)/(9*theta[d+1]^4),
#             10*(3*sqrt(5)*theta[d+1] + 5*(w1+w2))/(9*theta[d+1]^4),
#             25/(9*theta[d+1]^4))
#     
#     
#     return(exp((10*s + sqrt(5)*theta[d+1]*(w1 + w2 - 2*m))/theta[d+1]^2) *
#              (e3 %*% lambda31 * pnorm((muc - w1)/sqrt(s)) +
#                 e3 %*% lambda32 * sqrt(s)/sqrt(2*pi) * exp(-(w1 - muc)^2/(2*s))) +
#              exp(-(sqrt(5)*(w1-w2)/theta[d+1])) *
#              (e4 %*% lambda41 * (pnorm((w1 - m)/sqrt(s)) - pnorm((w2 - m)/sqrt(s))) +
#                 e4 %*% lambda42 * sqrt(s)/sqrt(2*pi) * exp(-(w2 - m)^2/(2*s)) -
#                 e4 %*% lambda43 * sqrt(s)/sqrt(2*pi) * exp(-(w1 - m)^2/(2*s)) ) +
#              exp((10*s - sqrt(5)*theta[d+1]*(w1 + w2 - 2*m))/theta[d+1]^2) *
#              (e5 %*% lambda51 * pnorm((w2 - mud)/sqrt(s)) +
#                 e5 %*% lambda52 * sqrt(s)/sqrt(2*pi) * exp(-(w2 - mud)^2/(2*s))))
#   }}