### Function to calculate IMSPE of n+1 when the medium fidelity is added ###
library(lhs)
library(laGP)
source("GP.R")
source("closed.R")

IMSPE2 <- function(x, newx, fit1, fit2, fit3, mc.sample=10000){ 
  # x; usually grid
  # newx; new point which will be added
  # fit1; first layer's emulator
  # fit2; second layer's emulator f_M(X2, f1(x2))
  # fit3; third layer's emulator f_H(X3, f_M(X3, f1(X3)))
  
  x1.sample <- rnorm(mc.sample, mean=pred.GP(fit1, newx)$mu, sd=sqrt(pred.GP(fit1, newx)$sig2))
  x2.sample <- rnorm(mc.sample, 
                     mean=pred.GP(fit2, cbind(newx, mean(x1.sample)))$mu, 
                     sd=sqrt(pred.GP(fit2, cbind(newx, mean(x1.sample)))$sig2))
  mu.cand <- mean(x2.sample) # f2(newx)
  
  fit2new <- GP(rbind(t(t(fit2$X)*attr(fit2$X, "scaled:scale")+attr(fit2$X, "scaled:center")), cbind(newx, mean(x1.sample))),
                rbind(fit2$y*attr(fit2$y, "scaled:scale")+attr(fit2$y, "scaled:center"), mu.cand))

  return(mean(closed2(x, fit1, fit2new, fit3)$sig2))
}
