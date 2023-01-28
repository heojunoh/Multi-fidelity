library(lhs)
library(laGP)
source("GP.R")
source("KOH.R")
source("closed.R")

### synthetic function ###
fl <- function(x, l){
  term1 <- sin(2*pi*x)
  term2 <- 0.2 * sin(8*pi*x)
  
  term1 + term2*0.8^l*5 # term1 + error term
}

### training data ###
n1 <- 30; n2 <- 20; n3 <- 8
set.seed(1)
X3 <- maximinLHS(n3, 1) # x^H
y3 <- fl(X3, l=5)
X2 <- matrix(c(X3, maximinLHS(n2-n3, 1))) # x^M
y2 <- fl(X2, l=3)
X1 <- matrix(c(X2, maximinLHS(n1-n2, 1))) # x^L
y1 <- fl(X1, l=1)

### model fitting for f1 ###
eps <- sqrt(.Machine$double.eps)
fit.GP1 <- GP(X1, y1)

### model fitting using (x2, f1(x2)) ###
w1.x2 <- pred.GP(fit.GP1, X2)$mu # can interpolate; nested
X2new <- cbind(X2, w1.x2) # combine (X2, f1(x2))
fit.GP2new <- GP(X2new, y2) # model fitting for f_M(X2, f1(x2))

### model fitting using (x3, f2(x3, f1(x3))) ###
w1.x3 <- pred.GP(fit.GP1, X3)$mu # can interpolate; nested
w2.x3 <- pred.GP(fit.GP2new, cbind(X3, w1.x3))$mu # can interpolate; nested
X3new <- cbind(X3, w2.x3) # combine (X3, f2(x3, f1(x3)))
fit.GP3new <- GP(X3new, y3) # model fitting for f_H(X3, f2(x3, f1(x3)))


### test data ###
x <- seq(0,1,0.001)


### closed ###
predy <- closed2(x, fit.GP1, fit.GP2new, fit.GP3new)$mu
predsig2 <- closed2(x, fit.GP1, fit.GP2new, fit.GP3new)$sig2


### compared to single fidelity ###
fit.GP3 <- GP(X3, y3)
pred3 <- pred.GP(fit.GP3, x)

### direct fitting; not using closed form. f1(u) and f_M(u) from (u, f_M(u, f1(u))) are random variables.
w1.x <- rnorm(length(x), mean=pred.GP(fit.GP1, x)$mu, sd=sqrt(pred.GP(fit.GP1, x)$sig2)) # sample f1(x)
w2.x <- rnorm(length(x), mean=pred.GP(fit.GP2new, cbind(x, w1.x))$mu, sd=sqrt(pred.GP(fit.GP2new, cbind(x, w1.x))$sig2)) # It makes wiggly.
xnew <- cbind(x, w2.x)
pred3new <- pred.GP(fit.GP3new, xnew) # not closed form

### RMSE ###
sqrt(mean((predy-fl(x, l=5))^2)) # closed form
sqrt(mean((pred3new$mu-fl(x, l=5))^2)) # not closed form
sqrt(mean((pred3$mu-fl(x, l=5))^2)) # single fidelity



#############
### IMSPE ###
#############
g <- seq(0,1,0.001) # more than 1-dim, use expand.grid()
mean(closed2(g, fit.GP1, fit.GP2new, fit.GP3new)$sig2) # current IMSPE

### Add 1 points to the low-fidelity data.
Icand1 <- c(rep(0, length(g))) # IMSPE candidates
Icand2 <- c(rep(0, length(g))) # IMSPE candidates
Icand3 <- c(rep(0, length(g))) # IMSPE candidates

for(i in 1:length(Icand1)){ # no true, no need to fit just pred
  Icand1[i] <- IMSPE1(g, g[i], fit.GP1, fit.GP2new, fit.GP3new)
}
for(i in 1:length(Icand2)){ # no true, no need to fit just pred
  Icand2[i] <- IMSPE2(g, g[i], fit.GP1, fit.GP2new, fit.GP3new)
}
for(i in 1:length(Icand3)){ # no true, no need to fit just pred
  Icand3[i] <- IMSPE3(g, g[i], fit.GP1, fit.GP2new, fit.GP3new)
}

plot(g, Icand1, type="l", lwd=2, col=3, ylim=range(Icand1))
plot(g, Icand2, type="l", lwd=2, col=3, ylim=range(Icand2))
plot(g, Icand3, type="l", lwd=2, col=3, ylim=range(Icand3))



# Idiff <- I current - low, medium, high



# Idiff * cost (C_L, C_L+C_M, C_L+C_M+C_H)








