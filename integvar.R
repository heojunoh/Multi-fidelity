
integvar <- function(x, fit, mc.sample=100){
  
  intvar1 <- c(rep(0, dim(t(x))[2])) # IMSPE candidates
  intvar2 <- c(rep(0, dim(t(x))[2])) # IMSPE candidates
  pseudointvar1 <- c(rep(0, mc.sample))
  pseudointvar2 <- c(rep(0, mc.sample))
  
  fit1 <- f1 <- fit$fit1
  fit2 <- f2 <- fit$fit2
  constant <- fit$constant
  kernel <- fit$kernel
  x <- matrix(x, nrow=dim(t(x))[2])
  
  x.center1 <- attr(fit1$X, "scaled:center")
  x.scale1 <- attr(fit1$X, "scaled:scale")
  y.center1 <- attr(fit1$y, "scaled:center")
  
  x.center2 <- attr(fit2$X, "scaled:center")
  x.scale2 <- attr(fit2$X, "scaled:scale")
  y.center2 <- attr(fit2$y, "scaled:center")
  
  
  for(i in 1:length(intvar1)){
    ### Choose level 1 ###
    newx <- matrix(x[i,], nrow=1)
    
    if(kernel=="sqex"){
      x1.sample <- rnorm(mc.sample, mean=pred.GP(f1, newx)$mu, sd=sqrt(pred.GP(f1, newx)$sig2))
    }else if(kernel=="matern1.5"){
      x1.sample <- rnorm(mc.sample, mean=pred.matGP(f1, newx)$mu, sd=sqrt(pred.matGP(f1, newx)$sig2))
    }else if(kernel=="matern2.5"){
      x1.sample <- rnorm(mc.sample, mean=pred.matGP(f1, newx)$mu, sd=sqrt(pred.matGP(f1, newx)$sig2))
    }
    
    newx1 <- matrix((newx-x.center1)/x.scale1, nrow=1)
    
    
    for(j in 1:mc.sample){
      
      if(kernel=="sqex"){
        x2.sample <- pred.GP(fit2, cbind(newx, x1.sample[j]))$mu
      }else if(kernel=="matern1.5"){
        x2.sample <- pred.matGP(fit2, cbind(newx, x1.sample[j]))$mu
      }else if(kernel=="matern2.5"){
        x2.sample <- pred.matGP(fit2, cbind(newx, x1.sample[j]))$mu
      }
      
      ### update Ki1
      if(kernel=="sqex"){
        v.next1 <- drop(covar.sep(X1=newx1, d=f1$theta, g=0) -
                          t(covar.sep(X1=f1$X, X2=newx1, d=f1$theta, g=0)) %*%
                          f1$Ki %*%
                          covar.sep(X1=f1$X, X2=newx1, d=f1$theta, g=0))
        g.next1 <- - drop(solve(v.next1)) * f1$Ki %*% covar.sep(X1=f1$X, X2=newx1, d=f1$theta, g=0)
      }else if(kernel=="matern1.5"){
        v.next1 <- drop(cor.sep(X=newx1, theta=f1$theta, nu=1.5) -
                          t(cor.sep(X=f1$X, x=newx1, theta=f1$theta, nu=1.5)) %*%
                          f1$Ki %*%
                          cor.sep(X=f1$X, x=newx1, theta=f1$theta, nu=1.5))
        g.next1 <- - drop(solve(v.next1)) * f1$Ki %*% cor.sep(X=f1$X, x=newx1, theta=f1$theta, nu=1.5)
      }else if(kernel=="matern2.5"){
        v.next1 <- drop(cor.sep(X=newx1, theta=f1$theta, nu=2.5) -
                          t(cor.sep(X=f1$X, x=newx1, theta=f1$theta, nu=2.5)) %*%
                          f1$Ki %*%
                          cor.sep(X=f1$X, x=newx1, theta=f1$theta, nu=2.5))
        g.next1 <- - drop(solve(v.next1)) * f1$Ki %*% cor.sep(X=f1$X, x=newx1, theta=f1$theta, nu=2.5)
      }
      
      fit1$Ki <- rbind(cbind(f1$Ki+g.next1%*%t(g.next1)*v.next1, g.next1),
                       cbind(t(g.next1), solve(v.next1)))
      
      fit1$X <- rbind(f1$X, newx1)
      attr(fit1$X, "scaled:center") <- x.center1
      attr(fit1$X, "scaled:scale") <- x.scale1
      
      if(constant){
        fit1$y <- c(f1$y, x1.sample[j])
      }else{
        fit1$y <- c(f1$y, x1.sample[j]-y.center1)
        attr(fit1$y, "scaled:center") <- y.center1
      }
      
      fit1$tau2hat <- drop(t(fit1$y) %*% fit1$Ki %*% fit1$y / length(fit1$y))
      
      fit$fit1 <- fit1
      
      pseudointvar1[j] <- mean(predRNAmf(fit, x)$sig2)
      
      ### Choose level 2 ###
      ### update Ki2
      newx2 <- t((t(cbind(newx, x1.sample[j]))-x.center2)/x.scale2)
      
      if(kernel=="sqex"){
        v.next2 <- drop(covar.sep(X1=newx2, d=f2$theta, g=0) -
                          t(covar.sep(X1=f2$X, X2=newx2, d=f2$theta, g=0)) %*%
                          f2$Ki %*%
                          covar.sep(X1=f2$X, X2=newx2, d=f2$theta, g=0))
        g.next2 <- - drop(solve(v.next2)) * f2$Ki %*% covar.sep(X1=f2$X, X2=newx2, d=f2$theta, g=0)
      }else if(kernel=="matern1.5"){
        v.next2 <- drop(cor.sep(X=newx2, theta=f2$theta, nu=1.5) -
                          t(cor.sep(X=f2$X, x=newx2, theta=f2$theta, nu=1.5)) %*%
                          f2$Ki %*%
                          cor.sep(X=f2$X, x=newx2, theta=f2$theta, nu=1.5))
        g.next2 <- - drop(solve(v.next2)) * f2$Ki %*% cor.sep(X=f2$X, x=newx2, theta=f2$theta, nu=1.5)
      }else if(kernel=="matern2.5"){
        v.next2 <- drop(cor.sep(X=newx2, theta=f2$theta, nu=2.5) -
                          t(cor.sep(X=f2$X, x=newx2, theta=f2$theta, nu=2.5)) %*%
                          f2$Ki %*%
                          cor.sep(X=f2$X, x=newx2, theta=f2$theta, nu=2.5))
        g.next2 <- - drop(solve(v.next2)) * f2$Ki %*% cor.sep(X=f2$X, x=newx2, theta=f2$theta, nu=2.5)
      }
      
      fit2$Ki <- rbind(cbind(f2$Ki+g.next2%*%t(g.next2)*v.next2, g.next2),
                       cbind(t(g.next2), solve(v.next2)))
      
      fit2$X <- rbind(f2$X, newx2)
      attr(fit2$X, "scaled:center") <- x.center2
      attr(fit2$X, "scaled:scale") <- x.scale2
      
      if(constant){
        fit2$y <- c(f2$y, x2.sample)
      }else{
        fit2$y <- c(f2$y, x2.sample-y.center2)
        attr(fit2$y, "scaled:center") <- y.center2
      }
      
      fit2$tau2hat <- drop(t(fit2$y) %*% fit2$Ki %*% fit2$y / length(fit2$y))
      
      fit$fit2 <- fit2
      
      pseudointvar2[j] <- mean(predRNAmf(fit, x)$sig2)
      
      fit$fit1 <- fit1 <- f1
      fit$fit2 <- fit2 <- f2
    }
    
    intvar1[i] <- mean(pseudointvar1)
    intvar2[i] <- mean(pseudointvar2)
  }
  
  return(list(intvar1=intvar1, intvar2=intvar2))
}
