### Function to calculate IMSPE of n+1 when the high fidelity is added ###
library(lhs)
library(laGP)
source("GP.R")
source("closed.R")

IMSPE3 <- function(x, newx, fit1, fit2, fit3, mc.sample=10000){ 
  # x; usually grid
  # newx; new point which will be added
  # fit1; first layer's emulator
  # fit2; second layer's emulator f_M(X2, f1(x2))
  # fit3; third layer's emulator f_H(X3, f_M(X3, f1(X3)))
  
  x1.sample <- rnorm(mc.sample, mean=pred.GP(fit1, newx)$mu, sd=sqrt(pred.GP(fit1, newx)$sig2))
  x2.sample <- rnorm(mc.sample, 
                     mean=pred.GP(fit2, cbind(newx, mean(x1.sample)))$mu, 
                     sd=sqrt(pred.GP(fit2, cbind(newx, mean(x1.sample)))$sig2))
  x3.sample <- rnorm(mc.sample, 
                     mean=pred.GP(fit3, cbind(newx, mean(x2.sample)))$mu, 
                     sd=sqrt(pred.GP(fit3, cbind(newx, mean(x2.sample)))$sig2))
  mu.cand <- mean(x3.sample) # f3(newx)
  
  fit3new <- GP(rbind(t(t(fit3$X)*attr(fit3$X, "scaled:scale")+attr(fit3$X, "scaled:center")), cbind(newx, mean(x2.sample))),
                rbind(fit3$y*attr(fit3$y, "scaled:scale")+attr(fit3$y, "scaled:center"), mu.cand))
  
  return(mean(closed2(x, fit1, fit2, fit3new)$sig2))
}
