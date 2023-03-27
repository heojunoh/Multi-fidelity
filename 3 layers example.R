library(lhs)
library(laGP)
source("GP.R")
source("KOH.R")
source("closed.R")
source("score.R")

### synthetic function ###
fl <- function(x, l){
  term1 <- sin(2*pi*x)
  term2 <- 0.2 * sin(8*pi*x)

  (term1 + term2*5 + (term1+term2)^3)*(1+0.8^l) 
}

### training data ###
n1 <- 11; n2 <- 9; n3 <- 7
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
x <- seq(0,1,0.01)


### closed ###
predy <- closed2(x, fit.GP1, fit.GP2new, fit.GP3new)$mu
predsig2 <- closed2(x, fit.GP1, fit.GP2new, fit.GP3new)$sig2


### KOH method ###
fit.KOH3 <- fit.KOH(X1, X2, X3, y1, y2, y3)
pred.KOH3 <- pred.KOH(fit.KOH3, x)
mx2 <- pred.KOH3$mu
koh.var2 <- pred.KOH3$sig2


### direct fitting; not using closed form. f1(u) and f_M(u) from (u, f_M(u, f1(u))) are random variables.
# w1.x <- rnorm(length(x), mean=pred.GP(fit.GP1, x)$mu, sd=sqrt(pred.GP(fit.GP1, x)$sig2)) # sample f1(x)
# w2.x <- rnorm(length(x), mean=pred.GP(fit.GP2new, cbind(x, w1.x))$mu, sd=sqrt(pred.GP(fit.GP2new, cbind(x, w1.x))$sig2))
# xxnew <- cbind(x, w2.x)
# pred3new <- pred.GP(fit.GP3new, xxnew) # not closed form
w1.x <- c(rep(NA, length(x)))
w2.x <- c(rep(NA, length(x)))
for(j in 1:length(x)){
  w1.x[j] <- mean(rnorm(10000, mean=pred.GP(fit.GP1, x[j])$mu, sd=sqrt(pred.GP(fit.GP1, x[j])$sig2)))
}
for(j in 1:length(x)){
  w2.x[j] <- mean(rnorm(10000, mean=pred.GP(fit.GP2new, cbind(x[j], w1.x[j]))$mu, sd=sqrt(pred.GP(fit.GP2new, cbind(x[j], w1.x[j]))$sig2)))
}

xxnew <- cbind(x, w2.x)
pred3new <- pred.GP(fit.GP3new, xxnew) # not closed form

### prediction of original GP with single fidelity ###
fit.GP3 <- GP(X3, y3)
pred3 <- pred.GP(fit.GP3, x)

### plot ###
plot(x, predy, type="l", lwd=2, col=3, 
     ylim=range(c(predy+1.96*sqrt(predsig2*length(y3)/(length(y3)-2)), predy-1.96*sqrt(predsig2*length(y3)/(length(y3)-2)), 
                  pred3new$mu+1.96*sqrt(pred3new$sig2*length(y3)/(length(y3)-2)), pred3new$mu-1.96*sqrt(pred3new$sig2*length(y3)/(length(y3)-2)),
                  mx2+1.96*sqrt(koh.var2*length(y3)/(length(y3)-2)), mx2-1.96*sqrt(koh.var2*length(y3)/(length(y3)-2))
     ))) # Green; Closed form
plot(x, predy, type="l", lwd=2, col=3, 
     ylim=range(c(predy+1.96*sqrt(predsig2*length(y3)/(length(y3)-2)), predy-1.96*sqrt(predsig2*length(y3)/(length(y3)-2))
     ))) 
lines(x, predy+1.96*sqrt(predsig2*length(y3)/(length(y3)-2)), col=3, lty=2)
lines(x, predy-1.96*sqrt(predsig2*length(y3)/(length(y3)-2)), col=3, lty=2)


lines(x, pred3new$mu, lwd=2, col=6) # Purple; Direct fitting
lines(x, pred3new$mu+1.96*sqrt(pred3new$sig2*length(y3)/(length(y3)-2)), col=6, lty=2)
lines(x, pred3new$mu-1.96*sqrt(pred3new$sig2*length(y3)/(length(y3)-2)), col=6, lty=2)

# lines(x, pred3$mu, lwd=2, col=4) # Blue; Single fidelity
# lines(x, pred3$mu+1.96*sqrt(pred3$sig2*length(y3)/(length(y3)-2)), col=4, lty=2)
# lines(x, pred3$mu-1.96*sqrt(pred3$sig2*length(y3)/(length(y3)-2)), col=4, lty=2)

lines(x, mx2, lwd=2, col=7) # Yellow; KOH
lines(x, mx2+1.96*sqrt(koh.var2*length(y3)/(length(y3)-2)), col=7, lty=2)
lines(x, mx2-1.96*sqrt(koh.var2*length(y3)/(length(y3)-2)), col=7, lty=2)

curve(fl(x,l=5),add=TRUE, col=1,lwd=2,lty=2) # high fidelity(TRUE); Black

### RMSE ###
sqrt(mean((pred3$mu-fl(x, l=5))^2)) # single fidelity
sqrt(mean((predy-fl(x, l=5))^2)) # closed form
sqrt(mean((pred3new$mu-fl(x, l=5))^2)) # not closed form
sqrt(mean((mx2-fl(x, l=5))^2)) # KOH

mean(score(fl(x, l=5), pred3$mu, pred3$sig2)) # single fidelity
mean(score(fl(x, l=5), predy, predsig2)) # closed form
mean(score(fl(x, l=5), pred3new$mu, pred3new$sig2)) # not closed form
mean(score(fl(x, l=5), mx2, koh.var2)) # KOH

median(score(fl(x, l=5), pred3$mu, pred3$sig2)) # single fidelity
median(score(fl(x, l=5), predy, predsig2)) # closed form
median(score(fl(x, l=5), pred3new$mu, pred3new$sig2)) # not closed form
median(score(fl(x, l=5), mx2, koh.var2)) # KOH

mean(crps(fl(x, l=5), pred3$mu, pred3$sig2)) # single fidelity
mean(crps(fl(x, l=5), predy, predsig2)) # closed form
mean(crps(fl(x, l=5), pred3new$mu, pred3new$sig2)) # not closed form
mean(crps(fl(x, l=5), mx2, koh.var2)) # KOH

median(crps(fl(x, l=5), pred3$mu, pred3$sig2)) # single fidelity
median(crps(fl(x, l=5), predy, predsig2)) # closed form
median(crps(fl(x, l=5), pred3new$mu, pred3new$sig2)) # not closed form
median(crps(fl(x, l=5), mx2, koh.var2)) # KOH


sum(predsig2)
sum(pred3new$sig2)
sum(koh.var2)

