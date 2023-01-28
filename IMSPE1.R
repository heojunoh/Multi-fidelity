### Function to calculate IMSPE of n+1 when the low fidelity is added ###
library(lhs)
library(laGP)
source("GP.R")
source("closed.R")

IMSPE1 <- function(x, newx, fit1, fit2, fit3, mc.sample=10000){ 
  # x; usually grid
  # newx; new point which will be added
  # fit1; first layer's emulator
  # fit2; second layer's emulator f_M(X2, f1(x2))
  # fit3; third layer's emulator f_H(X3, f_M(X3, f1(X3)))
  
  x1.sample <- rnorm(mc.sample, mean=pred.GP(fit1, newx)$mu, sd=sqrt(pred.GP(fit1, newx)$sig2))
  mu.cand <- mean(x1.sample) # f1(newx)

  fit1new <- GP(rbind(t(t(fit1$X)*attr(fit1$X, "scaled:scale")+attr(fit1$X, "scaled:center")), newx),
                rbind(fit1$y*attr(fit1$y, "scaled:scale")+attr(fit1$y, "scaled:center"), mu.cand))
  
  return(mean(closed2(x, fit1new, fit2, fit3)$sig2))
}
