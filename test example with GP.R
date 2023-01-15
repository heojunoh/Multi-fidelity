library(lhs)
library(laGP)
source("GP.R")

### synthetic function ###
fl <- function(x, l){
  term1 <- sin(2*pi*x)
  term2 <- 0.2 * sin(8*pi*x)
  
  term1 + term2*0.8^l*5 # term1 + error term
}

# curve(fl(x,l=5)) # highest fidelity
# curve(fl(x,l=3), add=TRUE, col=2, lty=2) # intermidiate fidelity
# curve(fl(x,l=1), add=TRUE, col=2, lty=3) # lowest fidelity

### training data ###
n1 <- 30; n2 <- 15; n3 <- 8 # When n1 is less than 16, the new one performs bad(why?)
X1 <- maximinLHS(n1, 1) # x^L
y1 <- fl(X1, l=1)
X2 <- maximinLHS(n2, 1) # x^M
y2 <- fl(X2, l=3)
X3 <- maximinLHS(n3, 1) # x^H
y3 <- fl(X3, l=5)

### test data ###
x <- seq(0,1,0.01)

### model fitting for f1 ###
eps <- sqrt(.Machine$double.eps)
fit.GP1 <- GP(X1, y1)
pred1 <- pred.GP(fit.GP1, x)

# plot(X1, y1, ylim=c(-2,2), col=2)
# lines(x, pred1$mu, lwd=2,col=2)
# lines(x, pred1$mu+1.96*sqrt(pred1$sig2*length(y1)/(length(y1)-2)), col=2, lty=2)
# lines(x, pred1$mu-1.96*sqrt(pred1$sig2*length(y1)/(length(y1)-2)), col=2, lty=2)
# curve(fl(x,l=1), add=TRUE, lwd=2, col=1) # low fidelity(TRUE)

### model fitting using (x2, f1(x2)) ###
fx1 <- pred.GP(fit.GP1, matrix(X2,ncol=1)) # f1(x2)

mc.sample <- 1000
pred3 <- matrix(0,ncol=length(x),nrow=mc.sample)
for(i in 1:mc.sample){
  w1 <- rnorm(n2, mean=fx1$mu, sd=sqrt(fx1$sig2)) # samples from f1(x2)
  X2new <- cbind(X2, w1) # combine the training data, (f_L(x^M), x^M)
  fit.GP2new <- GP(X2new, y2)
  W1.x <- rnorm(length(x), mean=pred1$mu, sd=sqrt(pred1$sig2)) # test data
  pred2new <- pred.GP(fit.GP2new, cbind(x, W1.x))

  pred3[i,] <- rnorm(length(x), mean=pred2new$mu, sd=sqrt(pred2new$sig2))
}

plot(x,colMeans(pred3), type="l", lwd=2, col=3, ylim=range(pred3)) # Green
lines(x,apply(pred3,2,quantile,0.975), col=3, lty=2) # upper bound
lines(x,apply(pred3,2,quantile,0.025), col=3, lty=2) # lower bound
points(X2, y2) # sample points


### compared to single fidelity ###
fit.GP2 <- GP(X2, y2)
pred2 <- pred.GP(fit.GP2, x)

lines(x, pred2$mu, lwd=2, col=4) # Blue
lines(x, pred2$mu+1.96*sqrt(pred2$sig2*length(y2)/(length(y2)-2)), col=4, lty=2)
lines(x, pred2$mu-1.96*sqrt(pred2$sig2*length(y2)/(length(y2)-2)), col=4, lty=2)

curve(fl(x,l=3),add=TRUE, col=1,lwd=2,lty=2) # medium fidelity(TRUE); Black

print(sqrt(mean((colMeans(pred3) - fl(x,l=3))^2))) # 0.0007
print(sqrt(mean((pred2$mu - fl(x,l=3))^2))) # 0.018
