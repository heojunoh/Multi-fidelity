### Function to calculate IMSPE of n+1 when the low fidelity is added ###
library(lhs)
library(laGP)
source("GP.R")
source("closed.R")

# IMSPE1 <- function(x, newx, fit1, fit2, fit3, mc.sample=10000){
#   ### This is fitting the emulator adding (x and y from MC) ###
#   # x; usually grid
#   # newx; new point which will be added
#   # fit1; first layer's emulator
#   # fit2; second layer's emulator f_M(X2, f1(x2))
#   # fit3; third layer's emulator f_H(X3, f_M(X3, f1(X3)))
#   
#   x1.sample <- rnorm(mc.sample, mean=pred.GP(fit1, newx)$mu, sd=sqrt(pred.GP(fit1, newx)$sig2))
#   mu.cand1 <- mean(x1.sample) # f1(newx)
#   
#   fit1new <- GP(rbind(t(t(fit1$X)*attr(fit1$X, "scaled:scale")+attr(fit1$X, "scaled:center")), newx),
#                 rbind(fit1$y*attr(fit1$y, "scaled:scale")+attr(fit1$y, "scaled:center"), mu.cand1))
#   
#   return(mean(closed2(x, fit1new, fit2, fit3)$sig2), fit1new=fit1new)
# }


IMSPE1fast <- function(x, newx, fit1, fit2, fit3, mc.sample=10000){ 
  ### This is updating Ki using Ki_{n+1}'s closed form ###
  # x; usually grid
  # newx; new point which will be added, should be n=1
  # fit1; first layer's emulator
  # fit2; second layer's emulator f_M(X2, f1(x2))
  # fit3; third layer's emulator f_H(X3, f_M(X3, f1(X3)))
  
  x1.sample <- rnorm(mc.sample, mean=pred.GP(fit1, newx)$mu, sd=sqrt(pred.GP(fit1, newx)$sig2))
  mu.cand1 <- mean(x1.sample)-attr(fit1$y,"scaled:center") # f1(newx)
  
  x.center1 <- attr(fit1$X, "scaled:center")
  x.scale1 <- attr(fit1$X, "scaled:scale")
  y.center1 <- attr(fit1$y, "scaled:center")
  # y.scale1 <- attr(fit1$y, "scaled:scale")
  
  newx <- matrix((newx-attr(fit1$X,"scaled:center"))/attr(fit1$X,"scaled:scale"))
  
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
  
  fit1$y <- rbind(fit1$y, mu.cand1)
  attr(fit1$y, "scaled:center") <- y.center1
  # attr(fit1$y, "scaled:scale") <- y.scale1
  
  fit1$tau2hat <- drop(t(fit1$y) %*% fit1$Ki %*% fit1$y / length(fit1$y))
  
  return(list(IMSPE=mean(closed2(x, fit1, fit2, fit3)$sig2), fit1new=fit1))
}


IMSPE1select <- function(x, newx, fit1, fit2, fit3, mc.sample=10000){ 
  ### This is updating Ki using Ki_{n+1}'s closed form ###
  # x; usually grid
  # newx; new point which will be added, should be n=1
  # fit1; first layer's emulator
  # fit2; second layer's emulator f_M(X2, f1(x2))
  # fit3; third layer's emulator f_H(X3, f_M(X3, f1(X3)))
  
  y1.select <- fl(newx, l=1)-attr(fit1$y,"scaled:center")
  
  x.center1 <- attr(fit1$X, "scaled:center")
  x.scale1 <- attr(fit1$X, "scaled:scale")
  y.center1 <- attr(fit1$y, "scaled:center")
  # y.scale1 <- attr(fit1$y, "scaled:scale")
  
  newx <- matrix((newx-attr(fit1$X,"scaled:center"))/attr(fit1$X,"scaled:scale"))
  
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
  
  fit1$y <- rbind(fit1$y, y1.select)
  attr(fit1$y, "scaled:center") <- y.center1
  # attr(fit1$y, "scaled:scale") <- y.scale1
  
  fit1$tau2hat <- drop(t(fit1$y) %*% fit1$Ki %*% fit1$y / length(fit1$y))
  
  return(list(IMSPE=mean(closed2(x, fit1, fit2, fit3)$sig2), fit1new=fit1, fit2new=fit2, fit3new=fit3))
}

