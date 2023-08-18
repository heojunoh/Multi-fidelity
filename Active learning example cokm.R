library(lhs)
library(laGP)
source("GP.R")
source("cokm.R")
source("closed.R")

### synthetic function ###
fl <- function(x, l){
  term1 <- sin(2*pi*x)
  term2 <- 0.2 * sin(8*pi*x)
  
  term1 + term2*0.8^l*5 + (term1+term2)^3 + exp(-2*term1*term2) #* (term1+term2)^3 # term1 + error term + interaction term
}

### training data ###
n1 <- 9; n2 <- 7; n3 <- 5
set.seed(1)
X1 <- maximinLHS(n1, 1)
X2 <- maximinLHS(n2, 1)
X3 <- maximinLHS(n3, 1)

NestDesign <- NestedDesignBuild(design = list(X1,X2,X3))

X1 <- NestDesign$PX
X2 <- ExtractNestDesign(NestDesign,2)
X3 <- ExtractNestDesign(NestDesign,3)

y1 <- fl(X1, l=1)
y2 <- fl(X2, l=3)
y3 <- fl(X3, l=5)

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


### cokm ###
fit.muficokm <- MuFicokm(formula = list(~1,~1,~1), MuFidesign = NestDesign, covtype="gauss",
                         # coef.trend = list(0,c(0,0),c(0,0)),
                         lower=eps, upper=0.1,
                         response = list(y1,y2,y3), nlevel = 3)
pred.muficokm <- predict(fit.muficokm, x, "SK", cov.compute=TRUE)

### compared to single fidelity ###
fit.GP3 <- GP(X3, y3)
pred3 <- pred.GP(fit.GP3, x)

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

### RMSE ###
sqrt(mean((pred.muficokm$mean-fl(x, l=5))^2)) # cokm
sqrt(mean((pred3new$mu-fl(x, l=5))^2)) # not closed form
sqrt(mean((pred3$mu-fl(x, l=5))^2)) # single fidelity


#############
### IMSPE ###
#############
g <- seq(0,1,0.01) # more than 1-dim, use expand.grid()
Icurrent <- mean(pred.muficokm$sig2) # current mean of sig2
Icurrent

### Add 1 points to the low-fidelity data ###
Icand1cokm <- c(rep(0, length(g))) # IMSPE candidates
Icand2cokm <- c(rep(0, length(g))) # IMSPE candidates
Icand3cokm <- c(rep(0, length(g))) # IMSPE candidates

for(i in 1:length(Icand1cokm)){ # no true, no need to fit just pred
  IMSPE <- IMSPEcokm(g, g[i], fit.muficokm)
  Icand1cokm[i] <- IMSPE$k1
  Icand2cokm[i] <- IMSPE$k2
  Icand3cokm[i] <- IMSPE$k3
}

plot(g, Icand1cokm, type="l", lwd=2, col=3, ylim=range(Icand1cokm))
plot(g, Icand2cokm, type="l", lwd=2, col=3, ylim=range(Icand2cokm))
plot(g, Icand3cokm, type="l", lwd=2, col=3, ylim=range(Icand3cokm))

which.min(Icand1cokm)
which.min(Icand2cokm)
which.min(Icand3cokm)

### cokm update; Equation 6.6. in Surrogates ###
### ALC; How much can be improved. Equation 6.6. in Surrogates ###
alccokm <- c(Icand1cokm[which.min(Icand1cokm)], Icand2cokm[which.min(Icand2cokm)], Icand3cokm[which.min(Icand3cokm)] )
alccokm

which.min(alccokm*c(2,(2+4),(2+4+8)))
alccokm*c(2,(2+4),(2+4+8))

chosen <- matrix(0, ncol=2)
chosen[1,1] <- which.min(alccokm*c(2,(2+4),(2+4+8)))
chosen[1,2] <- which.min(cbind(Icand1cokm, Icand2cokm, Icand3cokm)[,chosen[1,1]])


### Plotting the chosen point ###
plot(x, pred.muficokm$mean, type="l", lwd=2, col=3, 
     ylim=range(c(pred.muficokm$mean+1.96*sqrt(pred.muficokm$sig2*length(y3)/(length(y3)-2)), pred.muficokm$mean-1.96*sqrt(pred.muficokm$sig2*length(y3)/(length(y3)-2))
     ))) 
lines(x, pred.muficokm$mean+1.96*sqrt(pred.muficokm$sig2*length(y3)/(length(y3)-2)), col=3, lty=2)
lines(x, pred.muficokm$mean-1.96*sqrt(pred.muficokm$sig2*length(y3)/(length(y3)-2)), col=3, lty=2)

curve(fl(x,l=5),add=TRUE, col=1,lwd=2,lty=2) # high fidelity(TRUE); Black

points(X1, y1, pch="1", col="red")
points(X2, y2, pch="2", col="red")
points(X3, y3, pch="3", col="red")

text(g[which.min(Icand1cokm)], pred.muficokm$mean[which.min(Icand1cokm)], expression("1*"), col="red")
text(g[which.min(Icand2cokm)], pred.muficokm$mean[which.min(Icand2cokm)], expression("2*"), col="red")
text(g[which.min(Icand3cokm)], pred.muficokm$mean[which.min(Icand3cokm)], expression("3*"), col="red")



#################
### Add point ###
#################
Iselect <- IMSPEcokmselect(g, g[which.min(Icand1cokm)], fit.muficokm, level=1)

### closed ###
predy <- predict(Iselect$fit, x, "SK")$mean
predsig2 <- predict(Iselect$fit, x, "SK")$sig2

### RMSE ###  0.07955958 -> 0.07413344 -> 0.07678651 -> 
sqrt(mean((predy-fl(x, l=5))^2)) # closed form

#############
### IMSPE ###
#############
g <- seq(0,1,0.01) # more than 1-dim, use expand.grid()
Icurrent <- mean(predsig2) # current mean of sig2
Icurrent 
### Add 1 points to the low-fidelity data ###
Icand1cokm <- c(rep(0, length(g))) # IMSPE candidates
Icand2cokm <- c(rep(0, length(g))) # IMSPE candidates
Icand3cokm <- c(rep(0, length(g))) # IMSPE candidates

for(i in 1:length(Icand1cokm)){ # no true, no need to fit just pred
  if(any(chosen[,2]==i)){Icand1cokm[i] <- Icand2cokm[i] <- Icand3cokm[i] <- 0}else{
    IMSPE <- IMSPEcokm(g, g[i], Iselect$fit)
    Icand1cokm[i] <- IMSPE$k1
    Icand2cokm[i] <- IMSPE$k2
    Icand3cokm[i] <- IMSPE$k3
  }
}
if(any(Icand1cokm==0)){Icand1cokm[which(Icand1cokm==0)] <-  max(Icand1cokm)}
if(any(Icand2cokm==0)){Icand2cokm[which(Icand2cokm==0)] <-  max(Icand2cokm)}
if(any(Icand3cokm==0)){Icand3cokm[which(Icand3cokm==0)] <-  max(Icand3cokm)}

plot(g, Icand1cokm, type="l", lwd=2, col=3, ylim=range(Icand1cokm))
plot(g, Icand2cokm, type="l", lwd=2, col=3, ylim=range(Icand2cokm))
plot(g, Icand3cokm, type="l", lwd=2, col=3, ylim=range(Icand3cokm))

which.min(Icand1cokm)
which.min(Icand2cokm)
which.min(Icand3cokm)

### ALC; How much can be reduced ###
alccokm <- c(Icand1cokm[which.min(Icand1cokm)], Icand2cokm[which.min(Icand2cokm)], Icand3cokm[which.min(Icand3cokm)] )
alccokm

which.min(alccokm*c(2,(2+4),(2+4+8)))
alccokm*c(2,(2+4),(2+4+8))

chosen <- rbind(chosen, c(which.min(alccokm*c(2,(2+4),(2+4+8))), which.min(cbind(Icand1cokm, Icand2cokm, Icand3cokm)[,which.min(alccokm*c(2,(2+4),(2+4+8)))])))


### Plotting the chosen point ###
plot(x, predy, type="l", lwd=2, col=3, 
     ylim=range(c(predy+1.96*sqrt(predsig2*length(y3)/(length(y3)-2)), predy-1.96*sqrt(predsig2*length(y3)/(length(y3)-2))
     ))) 
lines(x, predy+1.96*sqrt(predsig2*length(y3)/(length(y3)-2)), col=3, lty=2)
lines(x, predy-1.96*sqrt(predsig2*length(y3)/(length(y3)-2)), col=3, lty=2)

curve(fl(x,l=5),add=TRUE, col=1,lwd=2,lty=2) # high fidelity(TRUE); Black

points(Iselect$fit$cok[[1]]@X, Iselect$fit$cok[[1]]@y, pch="1", col="red")
points(Iselect$fit$cok[[2]]@X, Iselect$fit$cok[[2]]@y, pch="2", col="red")
points(Iselect$fit$cok[[3]]@X, Iselect$fit$cok[[3]]@y, pch="3", col="red")

text(g[which.min(Icand1cokm)], predy[which.min(Icand1cokm)], expression("1*"), col="red")
text(g[which.min(Icand2cokm)], predy[which.min(Icand2cokm)], expression("2*"), col="red")
text(g[which.min(Icand3cokm)], predy[which.min(Icand3cokm)], expression("3*"), col="red")



#################
### Add point after 1 at 1###
#################
Iselect <- IMSPEcokmselect(g, g[which.min(Icand1cokm)], Iselect$fit, level=1)

### closed ###
predy <- predict(Iselect$fit, x, "SK")$mean
predsig2 <- predict(Iselect$fit, x, "SK")$sig2

### RMSE ###  0.07955958 -> 0.07413344 -> 0.07420309 -> 0.07138561 -> 0.07247886
sqrt(mean((predy-fl(x, l=5))^2)) # closed form

#############
### IMSPE ###
#############
g <- seq(0,1,0.01) # more than 1-dim, use expand.grid()
Icurrent <- mean(predsig2) # current mean of sig2
Icurrent 

### Add 1 points to the low-fidelity data ###
Icand1cokm <- c(rep(0, length(g))) # IMSPE candidates
Icand2cokm <- c(rep(0, length(g))) # IMSPE candidates
Icand3cokm <- c(rep(0, length(g))) # IMSPE candidates

for(i in 1:length(Icand1cokm)){ # no true, no need to fit just pred
  if(any(chosen[,2]==i)){Icand1cokm[i] <- Icand2cokm[i] <- Icand3cokm[i] <- 0}else{
    IMSPE <- IMSPEcokm(g, g[i], Iselect$fit)
    Icand1cokm[i] <- IMSPE$k1
    Icand2cokm[i] <- IMSPE$k2
    Icand3cokm[i] <- IMSPE$k3
  }
}
if(any(Icand1cokm==0)){Icand1cokm[which(Icand1cokm==0)] <-  max(Icand1cokm)}
if(any(Icand2cokm==0)){Icand2cokm[which(Icand2cokm==0)] <-  max(Icand2cokm)}
if(any(Icand3cokm==0)){Icand3cokm[which(Icand3cokm==0)] <-  max(Icand3cokm)}

plot(g, Icand1cokm, type="l", lwd=2, col=3, ylim=range(Icand1cokm))
plot(g, Icand2cokm, type="l", lwd=2, col=3, ylim=range(Icand2cokm))
plot(g, Icand3cokm, type="l", lwd=2, col=3, ylim=range(Icand3cokm))

which.min(Icand1cokm)
which.min(Icand2cokm)
which.min(Icand3cokm)

### ALC; How much can be reduced ###
alccokm <- c(Icand1cokm[which.min(Icand1cokm)], Icand2cokm[which.min(Icand2cokm)], Icand3cokm[which.min(Icand3cokm)] )
alccokm

which.min(alccokm*c(2,(2+4),(2+4+8)))
alccokm*c(2,(2+4),(2+4+8))

chosen <- rbind(chosen, c(which.min(alccokm*c(2,(2+4),(2+4+8))), which.min(cbind(Icand1cokm, Icand2cokm, Icand3cokm)[,which.min(alccokm*c(2,(2+4),(2+4+8)))])))


### Plotting the chosen point ###
plot(x, predy, type="l", lwd=2, col=3, 
     ylim=range(c(predy+1.96*sqrt(predsig2*length(y3)/(length(y3)-2)), predy-1.96*sqrt(predsig2*length(y3)/(length(y3)-2))
     ))) 
lines(x, predy+1.96*sqrt(predsig2*length(y3)/(length(y3)-2)), col=3, lty=2)
lines(x, predy-1.96*sqrt(predsig2*length(y3)/(length(y3)-2)), col=3, lty=2)

curve(fl(x,l=5),add=TRUE, col=1,lwd=2,lty=2) # high fidelity(TRUE); Black

points(Iselect$fit$cok[[1]]@X, Iselect$fit$cok[[1]]@y, pch="1", col="red")
points(Iselect$fit$cok[[2]]@X, Iselect$fit$cok[[2]]@y, pch="2", col="red")
points(Iselect$fit$cok[[3]]@X, Iselect$fit$cok[[3]]@y, pch="3", col="red")

text(g[which.min(Icand1cokm)], predy[which.min(Icand1cokm)], expression("1*"), col="red")
text(g[which.min(Icand2cokm)], predy[which.min(Icand2cokm)], expression("2*"), col="red")
text(g[which.min(Icand3cokm)], predy[which.min(Icand3cokm)], expression("3*"), col="red")




#################
### Add point at level 2 ###
#################
Iselect <- IMSPEcokmselect(g, g[which.min(Icand2cokm)], Iselect$fit, level=2)

### closed ###
predy <- predict(Iselect$fit, x, "SK")$mean
predsig2 <- predict(Iselect$fit, x, "SK")$sig2

### RMSE ###  0.07955958 -> 0.07413344 -> 0.07420309 -> 0.07138561 -> 0.07247886
sqrt(mean((predy-fl(x, l=5))^2)) # closed form

#############
### IMSPE ###
#############
g <- seq(0,1,0.01) # more than 1-dim, use expand.grid()
Icurrent <- mean(predsig2) # current mean of sig2
Icurrent 

### Add 1 points to the low-fidelity data ###
Icand1cokm <- c(rep(0, length(g))) # IMSPE candidates
Icand2cokm <- c(rep(0, length(g))) # IMSPE candidates
Icand3cokm <- c(rep(0, length(g))) # IMSPE candidates

for(i in 1:length(Icand1cokm)){ # no true, no need to fit just pred
  if(any(chosen[,2]==i)){Icand1cokm[i] <- Icand2cokm[i] <- Icand3cokm[i] <- 0}else{
    IMSPE <- IMSPEcokm(g, g[i], Iselect$fit)
    Icand1cokm[i] <- IMSPE$k1
    Icand2cokm[i] <- IMSPE$k2
    Icand3cokm[i] <- IMSPE$k3
  }
}
if(any(Icand1cokm==0)){Icand1cokm[which(Icand1cokm==0)] <-  max(Icand1cokm)}
if(any(Icand2cokm==0)){Icand2cokm[which(Icand2cokm==0)] <-  max(Icand2cokm)}
if(any(Icand3cokm==0)){Icand3cokm[which(Icand3cokm==0)] <-  max(Icand3cokm)}

plot(g, Icand1cokm, type="l", lwd=2, col=3, ylim=range(Icand1cokm))
plot(g, Icand2cokm, type="l", lwd=2, col=3, ylim=range(Icand2cokm))
plot(g, Icand3cokm, type="l", lwd=2, col=3, ylim=range(Icand3cokm))

which.min(Icand1cokm)
which.min(Icand2cokm)
which.min(Icand3cokm)

### ALC; How much can be reduced ###
alccokm <- c(Icand1cokm[which.min(Icand1cokm)], Icand2cokm[which.min(Icand2cokm)], Icand3cokm[which.min(Icand3cokm)] )
alccokm

which.min(alccokm*c(2,(2+4),(2+4+8)))
alccokm*c(2,(2+4),(2+4+8))

chosen <- rbind(chosen, c(which.min(alccokm*c(2,(2+4),(2+4+8))), which.min(cbind(Icand1cokm, Icand2cokm, Icand3cokm)[,which.min(alccokm*c(2,(2+4),(2+4+8)))])))


### Plotting the chosen point ###
plot(x, predy, type="l", lwd=2, col=3, 
     ylim=range(c(predy+1.96*sqrt(predsig2*length(y3)/(length(y3)-2)), predy-1.96*sqrt(predsig2*length(y3)/(length(y3)-2))
     ))) 
lines(x, predy+1.96*sqrt(predsig2*length(y3)/(length(y3)-2)), col=3, lty=2)
lines(x, predy-1.96*sqrt(predsig2*length(y3)/(length(y3)-2)), col=3, lty=2)

curve(fl(x,l=5),add=TRUE, col=1,lwd=2,lty=2) # high fidelity(TRUE); Black

points(Iselect$fit$cok[[1]]@X, Iselect$fit$cok[[1]]@y, pch="1", col="red")
points(Iselect$fit$cok[[2]]@X, Iselect$fit$cok[[2]]@y, pch="2", col="red")
points(Iselect$fit$cok[[3]]@X, Iselect$fit$cok[[3]]@y, pch="3", col="red")

text(g[which.min(Icand1cokm)], predy[which.min(Icand1cokm)], expression("1*"), col="red")
text(g[which.min(Icand2cokm)], predy[which.min(Icand2cokm)], expression("2*"), col="red")
text(g[which.min(Icand3cokm)], predy[which.min(Icand3cokm)], expression("3*"), col="red")





#################
### Add point at level 3 ###
#################
Iselect <- IMSPEcokmselect(g, g[which.min(Icand3cokm)], Iselect$fit, level=3)

### closed ###
predy <- predict(Iselect$fit, x, "SK")$mean
predsig2 <- predict(Iselect$fit, x, "SK")$sig2

### RMSE ###  0.07955958 -> 0.07413344 -> 0.07420309 -> 0.07138561 -> 0.07247886
sqrt(mean((predy-fl(x, l=5))^2)) # closed form

#############
### IMSPE ###
#############
g <- seq(0,1,0.01) # more than 1-dim, use expand.grid()
Icurrent <- mean(predsig2) # current mean of sig2
Icurrent 

### Add 1 points to the low-fidelity data ###
Icand1cokm <- c(rep(0, length(g))) # IMSPE candidates
Icand2cokm <- c(rep(0, length(g))) # IMSPE candidates
Icand3cokm <- c(rep(0, length(g))) # IMSPE candidates

for(i in 1:length(Icand1cokm)){ # no true, no need to fit just pred
  if(any(chosen[,2]==i)){Icand1cokm[i] <- Icand2cokm[i] <- Icand3cokm[i] <- 0}else{
    IMSPE <- IMSPEcokm(g, g[i], Iselect$fit)
    Icand1cokm[i] <- IMSPE$k1
    Icand2cokm[i] <- IMSPE$k2
    Icand3cokm[i] <- IMSPE$k3
  }
}
if(any(Icand1cokm==0)){Icand1cokm[which(Icand1cokm==0)] <-  max(Icand1cokm)}
if(any(Icand2cokm==0)){Icand2cokm[which(Icand2cokm==0)] <-  max(Icand2cokm)}
if(any(Icand3cokm==0)){Icand3cokm[which(Icand3cokm==0)] <-  max(Icand3cokm)}

plot(g, Icand1cokm, type="l", lwd=2, col=3, ylim=range(Icand1cokm))
plot(g, Icand2cokm, type="l", lwd=2, col=3, ylim=range(Icand2cokm))
plot(g, Icand3cokm, type="l", lwd=2, col=3, ylim=range(Icand3cokm))

which.min(Icand1cokm)
which.min(Icand2cokm)
which.min(Icand3cokm)

### ALC; How much can be reduced ###
alccokm <- c(Icand1cokm[which.min(Icand1cokm)], Icand2cokm[which.min(Icand2cokm)], Icand3cokm[which.min(Icand3cokm)] )
alccokm

which.min(alccokm*c(2,(2+4),(2+4+8)))
alccokm*c(2,(2+4),(2+4+8))

chosen <- rbind(chosen, c(which.min(alccokm*c(2,(2+4),(2+4+8))), which.min(cbind(Icand1cokm, Icand2cokm, Icand3cokm)[,which.min(alccokm*c(2,(2+4),(2+4+8)))])))


### Plotting the chosen point ###
plot(x, predy, type="l", lwd=2, col=3, 
     ylim=range(c(predy+1.96*sqrt(predsig2*length(y3)/(length(y3)-2)), predy-1.96*sqrt(predsig2*length(y3)/(length(y3)-2))
     ))) 
lines(x, predy+1.96*sqrt(predsig2*length(y3)/(length(y3)-2)), col=3, lty=2)
lines(x, predy-1.96*sqrt(predsig2*length(y3)/(length(y3)-2)), col=3, lty=2)

curve(fl(x,l=5),add=TRUE, col=1,lwd=2,lty=2) # high fidelity(TRUE); Black

points(Iselect$fit$cok[[1]]@X, Iselect$fit$cok[[1]]@y, pch="1", col="red")
points(Iselect$fit$cok[[2]]@X, Iselect$fit$cok[[2]]@y, pch="2", col="red")
points(Iselect$fit$cok[[3]]@X, Iselect$fit$cok[[3]]@y, pch="3", col="red")

text(g[which.min(Icand1cokm)], predy[which.min(Icand1cokm)], expression("1*"), col="red")
text(g[which.min(Icand2cokm)], predy[which.min(Icand2cokm)], expression("2*"), col="red")
text(g[which.min(Icand3cokm)], predy[which.min(Icand3cokm)], expression("3*"), col="red")



