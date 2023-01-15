library(lhs)
library(laGP)
install.packages("mlegp")
library(mlegp)

### synthetic function ###
fl <- function(x, l){
  term1 <- sin(2*pi*x)
  term2 <- 0.2 * sin(8*pi*x)
  
  term1 + term2*0.8^l*5 # term1 + error term
}

### training data ###
n1 <- 30; n2 <- 15; 
X1 <- maximinLHS(n1, 1) # x^L
y1 <- fl(X1, l=1)
X2 <- maximinLHS(n2, 1) # x^M
y2 <- fl(X2, l=3)

### model fitting for f1 ###
eps <- sqrt(.Machine$double.eps)
d1 <- darg(list(mle=TRUE), X1)
gpi1 <- newGP(X1, y1, d=d1$start, g=eps, dK=TRUE)
mle1 <- mleGP(gpi1, param = "d", d1$min, d1$max)
pred1 <- predGP(gpi1, matrix(x,ncol=1), lite=TRUE)

### model fitting using (x2, f1(x2)) ###
X2new <- cbind(X2, predGP(gpi1, X2)$mean) # combine (X2, f1(x2))
d2new <- darg(list(mle=TRUE), X2new)
gpi2new <- newGPsep(X2new, y2, d=rep(d2new$start,2), g=eps, dK=TRUE)
mle2new <- mleGPsep(gpi2new, param = "d", rep(d2new$min,2), rep(d2new$max,2))

fit.GP2new <- createGP(X2new, y2, 1/mle2new$d, rep(2,ncol(X2new)), meanReg=1, sig2=mlegp(X2new, y2)$sig2, nugget=eps)


### test data ###
x <- seq(0,1,0.01)
xnew <- cbind(x, pred1$mean)

### calculate the closed form ###
Ci <- fit.GP2new$invVarMatrix
a <- Ci %*% y2
mu <- predGPsep(gpi2new, matrix(xnew, ncol=2), lite=TRUE)$mean
sig2 <- diag(predGPsep(gpi2new, xnew)$Sigma)
tau2hat <- drop(t(y2) %*% Ci %*% y2 / length(y2))

# mean
predy <- c(rep(0, length(x)))

for(j in 1: length(x)){
  predv <- c(rep(0, length(y2)))
  for(i in 1: length(y2)){
    predv[i] <- exp(-(xnew[j,1]-X2new[i,1])^2/mle2new$d[1]) *
      1/sqrt(1+2*sig2[j]/mle2new$d[2]) *
      exp(-(mu[j]-X2new[i,2])^2/(mle2new$d[2]+2*sig2[j]))
  }
  predy[j] <- sum(predv*a)
}

# var
predsig2 <- c(rep(0, length(x)))

for(k in 1: length(x)){
  mat1 <- matrix(0, length(y2), length(y2))
  for(i in 1: length(y2)){
    for(j in 1: length(y2)){
      mat1[i,j] <- exp(-((xnew[k,1]-X2new[i,1])/mle2new$d[1])^2-((xnew[k,1]-X2new[j,1])/mle2new$d[1])^2) 
    }}
  
  mat2 <- matrix(0, length(y2), length(y2))
  for(i in 1: length(y2)){
    for(j in 1: length(y2)){
      mat2[i,j] <- 1/sqrt(1+4*sig2[k]/mle2new$d[2]^2) * 
        exp(-((X2new[j,1]+X2new[j,1])/2-xnew[k,2])^2/(mle2new$d[2]^2/2+2*sig2[k])) * 
        exp(-(X2new[i,1]-X2new[j,1])^2/(2*mle2new$d[2]^2))
    }}
  
  predsig2[k] <- tau2hat + xnew[k,2]^2 - predy[k]^2 + sum((a %*% t(a) - tau2hat*Ci) * mat1 * mat2)
}

### compared to single fidelity ###
d2 <- darg(list(mle=TRUE), X2)
gpi2 <- newGP(X2, y2, d=d2$start, g=eps, dK=TRUE)
mle2 <- mleGP(gpi2, param = "d", d2$min, d2$max)
pred2 <- predGP(gpi2, matrix(x,ncol=1), lite=TRUE)

pred2new <- predGPsep(gpi2new, xnew) # not closed form

### delete gpi ###
deleteGP(gpi1)
deleteGPsep(gpi2new)
deleteGP(gpi2)

### RMSE ###
sqrt(mean((predy-fl(x, l=3))^2)) 
sqrt(mean((pred2$mean-fl(x, l=3))^2)) 
sqrt(mean((pred2new$mean-fl(x, l=3))^2)) 

### plot ###
plot(x,predy, type="l", lwd=2, col=3, ylim=c(-1.5,1.5)) # Green
lines(x, predy+1.96*sqrt(predsig2*length(y2)/(length(y2)-2)), col=3, lty=2)
lines(x, predy-1.96*sqrt(predsig2*length(y2)/(length(y2)-2)), col=3, lty=2)

lines(x, pred2$mean, lwd=2, col=4) # Blue
lines(x, pred2$mean+1.96*sqrt(pred2$s2*length(y2)/(length(y2)-2)), col=4, lty=2)
lines(x, pred2$mean-1.96*sqrt(pred2$s2*length(y2)/(length(y2)-2)), col=4, lty=2)

lines(x, pred2new$mean, lwd=2, col=7) # Yellow
lines(x, pred2new$mean+1.96*sqrt(diag(pred2new$Sigma)*length(y2)/(length(y2)-2)), col=7, lty=2)
lines(x, pred2new$mean-1.96*sqrt(diag(pred2new$Sigma)*length(y2)/(length(y2)-2)), col=7, lty=2)

curve(fl(x,l=3),add=TRUE, col=1,lwd=2,lty=2) # medium fidelity(TRUE); Black

