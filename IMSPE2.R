### Function to calculate IMSPE of n+1 when the medium fidelity is added ###
library(lhs)
library(laGP)
source("GP.R")
source("closed.R")

# IMSPE2 <- function(x, newx, fit1, fit2, fit3, mc.sample=10000){ 
#   # x; usually grid
#   # newx; new point which will be added
#   # fit1; first layer's emulator
#   # fit2; second layer's emulator f_M(X2, f1(x2))
#   # fit3; third layer's emulator f_H(X3, f_M(X3, f1(X3)))
#   x1.sample <- rnorm(mc.sample, mean=pred.GP(fit1, newx)$mu, sd=sqrt(pred.GP(fit1, newx)$sig2))
#   mu.cand1 <- mean(x1.sample) # f1(newx)
#   
#   fit1new <- GP(rbind(t(t(fit1$X)*attr(fit1$X, "scaled:scale")+attr(fit1$X, "scaled:center")), newx),
#                 rbind(fit1$y*attr(fit1$y, "scaled:scale")+attr(fit1$y, "scaled:center"), mu.cand1))
# 
#   x2.sample <- rnorm(mc.sample, 
#                      mean=pred.GP(fit2, cbind(newx, mu.cand1))$mu, 
#                      sd=sqrt(pred.GP(fit2, cbind(newx, mu.cand1))$sig2))
#   mu.cand2 <- mean(x2.sample) # f2(newx)
#   
#   fit2new <- GP(rbind(t(t(fit2$X)*attr(fit2$X, "scaled:scale")+attr(fit2$X, "scaled:center")), cbind(newx, mean(x1.sample))),
#                 rbind(fit2$y*attr(fit2$y, "scaled:scale")+attr(fit2$y, "scaled:center"), mu.cand2))
# 
#   return(mean(closed2(x, fit1new, fit2new, fit3)$sig2), fit1new=fit1new, fit2new=fit2new)
# }


IMSPE2fast <- function(x, newx, fit1, fit2, fit3, mc.sample=10000, constant=FALSE){ 
  ### This is updating Ki using Ki_{n+1}'s closed form ###
  # x; usually grid
  # newx; new point which will be added, should be n=1
  # fit1; first layer's emulator
  # fit2; second layer's emulator f_M(X2, f1(x2))
  # fit3; third layer's emulator f_H(X3, f_M(X3, f1(X3)))
  
  x1.sample <- rnorm(mc.sample, mean=pred.GP(fit1, newx)$mu, sd=sqrt(pred.GP(fit1, newx)$sig2))
  mu.cand1 <- mean(x1.sample) # f1(newx)
  
  x2.sample <- rnorm(mc.sample, 
                     mean=pred.GP(fit2, cbind(newx, mu.cand1))$mu, 
                     sd=sqrt(pred.GP(fit2, cbind(newx, mu.cand1))$sig2))
  mu.cand2 <- mean(x2.sample) # f2(newx)
  
  x.center1 <- attr(fit1$X, "scaled:center")
  x.scale1 <- attr(fit1$X, "scaled:scale")
  y.center1 <- attr(fit1$y, "scaled:center")
  # y.scale1 <- attr(fit1$y, "scaled:scale")
  
  x.center2 <- attr(fit2$X, "scaled:center")
  x.scale2 <- attr(fit2$X, "scaled:scale")
  y.center2 <- attr(fit2$y, "scaled:center")
  # y.scale2 <- attr(fit2$y, "scaled:scale")
  
  newx1 <- matrix((newx-attr(fit1$X,"scaled:center"))/attr(fit1$X,"scaled:scale")) 
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
    fit1$y <- rbind(fit1$y, mu.cand1)
  }else{
    fit1$y <- rbind(fit1$y, mu.cand1-attr(fit1$y,"scaled:center"))
    attr(fit1$y, "scaled:center") <- y.center1
    # attr(fit1$y, "scaled:scale") <- y.scale1
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
    fit2$y <- rbind(fit2$y, mu.cand2)
  }else{
    fit2$y <- rbind(fit2$y, mu.cand2-attr(fit2$y,"scaled:center"))
    attr(fit2$y, "scaled:center") <- y.center2
    # attr(fit2$y, "scaled:scale") <- y.scale2
  }
  
  fit2$tau2hat <- drop(t(fit2$y) %*% fit2$Ki %*% fit2$y / length(fit2$y))
  
  if(constant){
    return(list(IMSPE=mean(closed2(x, fit1, fit2, fit3, constant=TRUE)$sig2), fit1new=fit1, fit2new=fit2))
  }else{
    return(list(IMSPE=mean(closed2(x, fit1, fit2, fit3)$sig2), fit1new=fit1, fit2new=fit2))
  }
}


IMSPE2select <- function(x, newx, fit1, fit2, fit3, mc.sample=10000, constant=FALSE){ 
  ### This is updating Ki using Ki_{n+1}'s closed form ###
  # x; usually grid
  # newx; new point which will be added, should be n=1
  # fit1; first layer's emulator
  # fit2; second layer's emulator f_M(X2, f1(x2))
  # fit3; third layer's emulator f_H(X3, f_M(X3, f1(X3)))
  
  y1.select <- fl(newx, l=1)
  y2.select <- fl(newx, l=3)
  
  x.center1 <- attr(fit1$X, "scaled:center")
  x.scale1 <- attr(fit1$X, "scaled:scale")
  y.center1 <- attr(fit1$y, "scaled:center")
  # y.scale1 <- attr(fit1$y, "scaled:scale")
  
  x.center2 <- attr(fit2$X, "scaled:center")
  x.scale2 <- attr(fit2$X, "scaled:scale")
  y.center2 <- attr(fit2$y, "scaled:center")
  # y.scale2 <- attr(fit2$y, "scaled:scale")
  
  newx1 <- matrix((newx-attr(fit1$X,"scaled:center"))/attr(fit1$X,"scaled:scale")) 
  newx2 <- t((t(cbind(newx, y1.select))-attr(fit2$X,"scaled:center"))/attr(fit2$X,"scaled:scale"))
  
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
    fit1$y <- rbind(fit1$y, y1.select)
  }else{
    fit1$y <- rbind(fit1$y, y1.select-attr(fit1$y,"scaled:center"))
    attr(fit1$y, "scaled:center") <- y.center1
    # attr(fit1$y, "scaled:scale") <- y.scale1
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
    fit2$y <- rbind(fit2$y, y2.select)
  }else{
    fit2$y <- rbind(fit2$y, y2.select-attr(fit2$y,"scaled:center"))
    attr(fit2$y, "scaled:center") <- y.center2
    # attr(fit2$y, "scaled:scale") <- y.scale2
  }
  
  fit2$tau2hat <- drop(t(fit2$y) %*% fit2$Ki %*% fit2$y / length(fit2$y))

  if(constant){
    return(list(IMSPE=mean(closed2(x, fit1, fit2, fit3, constant=TRUE)$sig2), fit1new=fit1, fit2new=fit2, fit3new=fit3))
  }else{
    return(list(IMSPE=mean(closed2(x, fit1, fit2, fit3)$sig2), fit1new=fit1, fit2new=fit2, fit3new=fit3))
  }
}

