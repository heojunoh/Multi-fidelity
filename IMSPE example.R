library(lhs)
library(laGP)
source("GP.R")
source("KOH.R")
source("closed.R")
source("IMSPE1.R")
source("IMSPE2.R")
source("IMSPE3.R")

### synthetic function ###
fl <- function(x, l){
  term1 <- sin(2*pi*x)
  term2 <- 0.2 * sin(8*pi*x)
  
  term1 + term2*0.8^l*5 + (term1+term2)^3 + exp(-2*term1*term2) #* (term1+term2)^3 # term1 + error term + interaction term
}

### training data ###
n1 <- 9; n2 <- 7; n3 <- 5
set.seed(8)
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
fit.GP1 <- GP(X1, y1, constant=TRUE)

### model fitting using (x2, f1(x2)) ###
w1.x2 <- pred.GP(fit.GP1, X2)$mu # can interpolate; nested
X2new <- cbind(X2, w1.x2) # combine (X2, f1(x2))
fit.GP2new <- GP(X2new, y2, constant=TRUE) # model fitting for f_M(X2, f1(x2))

### model fitting using (x3, f2(x3, f1(x3))) ###
w1.x3 <- pred.GP(fit.GP1, X3)$mu # can interpolate; nested
w2.x3 <- pred.GP(fit.GP2new, cbind(X3, w1.x3))$mu # can interpolate; nested
X3new <- cbind(X3, w2.x3) # combine (X3, f2(x3, f1(x3)))
fit.GP3new <- GP(X3new, y3, constant=TRUE) # model fitting for f_H(X3, f2(x3, f1(x3)))


### test data ###
x <- seq(0,1,0.01)


### closed ###
predy <- closed2(x, fit.GP1, fit.GP2new, fit.GP3new, constant=TRUE)$mu
predsig2 <- closed2(x, fit.GP1, fit.GP2new, fit.GP3new, constant=TRUE)$sig2


### compared to single fidelity ###
fit.GP3 <- GP(X3, y3, constant=TRUE)
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
sqrt(mean((predy-fl(x, l=Inf))^2)) # closed form
sqrt(mean((pred3new$mu-fl(x, l=Inf))^2)) # not closed form
sqrt(mean((pred3$mu-fl(x, l=Inf))^2)) # single fidelity


#############
### IMSPE ###
#############
g <- seq(0,1,0.01) # more than 1-dim, use expand.grid()
Icurrent <- mean(closed2(g, fit.GP1, fit.GP2new, fit.GP3new, constant=TRUE)$sig2) # current IMSPE
Icurrent

### Add 1 points to the low-fidelity data ###
Icand1fast <- c(rep(0, length(g))) # IMSPE candidates
Icand2fast <- c(rep(0, length(g))) # IMSPE candidates
Icand3fast <- c(rep(0, length(g))) # IMSPE candidates

for(i in 1:length(Icand1fast)){ # no true, no need to fit just pred
  Icand1fast[i] <- IMSPE1fast(g, g[i], fit.GP1, fit.GP2new, fit.GP3new, constant=TRUE)$IMSPE
}
for(i in 1:length(Icand2fast)){ # no true, no need to fit just pred
  Icand2fast[i] <- IMSPE2fast(g, g[i], fit.GP1, fit.GP2new, fit.GP3new, constant=TRUE)$IMSPE
}
for(i in 1:length(Icand3fast)){ # no true, no need to fit just pred
  Icand3fast[i] <- IMSPE3fast(g, g[i], fit.GP1, fit.GP2new, fit.GP3new, constant=TRUE)$IMSPE
}

plot(g, Icand1fast, type="l", lwd=2, col=3, ylim=range(Icand1fast))
plot(g, Icand2fast, type="l", lwd=2, col=3, ylim=range(Icand2fast))
plot(g, Icand3fast, type="l", lwd=2, col=3, ylim=range(Icand3fast))

which.min(Icand1fast)
which.min(Icand2fast)
which.min(Icand3fast)

### Fast update; Equation 6.6. in Surrogates ###
### ALC; How much can be improved. Equation 6.6. in Surrogates ###
alcfast <- c(Icurrent - Icand1fast[which.min(Icand1fast)], Icurrent - Icand2fast[which.min(Icand2fast)], Icurrent - Icand3fast[which.min(Icand3fast)] )
alcfast

### cost; 1, 2, 3 ###
which.max(alcfast/c(2,(2+8),(2+8+32)))
alcfast/c(2,(2+8),(2+8+32))


chosen <- matrix(0, ncol=2)
chosen[1,1] <- which.max(alcfast/c(2,(2+8),(2+8+32)))
chosen[1,2] <- which.min(cbind(Icand1fast, Icand2fast, Icand3fast)[,chosen[1,1]])


### Plotting the chosen point ###
plot(x, predy, type="l", lwd=2, col=3, 
     ylim=range(c(predy+1.96*sqrt(predsig2*length(y3)/(length(y3)-2)), predy-1.96*sqrt(predsig2*length(y3)/(length(y3)-2))
     ))) 
lines(x, predy+1.96*sqrt(predsig2*length(y3)/(length(y3)-2)), col=3, lty=2)
lines(x, predy-1.96*sqrt(predsig2*length(y3)/(length(y3)-2)), col=3, lty=2)

curve(fl(x,l=5),add=TRUE, col=1,lwd=2,lty=2) # high fidelity(TRUE); Black

points(X1, y1, pch="1", col="red")
points(X2, y2, pch="2", col="red")
points(X3, y3, pch="3", col="red")

text(g[which.min(Icand1fast)], predy[which.min(Icand1fast)], expression("1*"), col="red")
text(g[which.min(Icand2fast)], predy[which.min(Icand2fast)], expression("2*"), col="red")
text(g[which.min(Icand3fast)], predy[which.min(Icand3fast)], expression("3*"), col="red")



#################
### Add point ###
#################
Iselect <- IMSPE1select(g, g[which.min(Icand1fast)], fit.GP1, fit.GP2new, fit.GP3new, constant=TRUE)

### closed ###
predy <- closed2(x, Iselect$fit1new, Iselect$fit2new, Iselect$fit3new, constant=TRUE)$mu
predsig2 <- closed2(x, Iselect$fit1new, Iselect$fit2new, Iselect$fit3new, constant=TRUE)$sig2

### RMSE ###  0.07955958 -> 0.07413344 -> 0.07678651 -> 
sqrt(mean((predy-fl(x, l=Inf))^2)) # closed form

#############
### IMSPE ###
#############
g <- seq(0,1,0.01) # more than 1-dim, use expand.grid()
Icurrent <- mean(closed2(g, Iselect$fit1new, Iselect$fit2new, Iselect$fit3new, constant=TRUE)$sig2) # current IMSPE, 0.002260122 -> 0.00084703 -> 0.001532291 -> 
Icurrent 
### Add 1 points to the low-fidelity data ###
Icand1fast <- c(rep(0, length(g))) # IMSPE candidates
Icand2fast <- c(rep(0, length(g))) # IMSPE candidates
Icand3fast <- c(rep(0, length(g))) # IMSPE candidates

for(i in 1:length(Icand1fast)){ # no true, no need to fit just pred
  if(any(chosen[,2]==i)){Icand1fast[i] <- 0}else{
    Icand1fast[i] <- IMSPE1fast(g, g[i], Iselect$fit1new, Iselect$fit2new, Iselect$fit3new, constant=TRUE)$IMSPE
  }
}
for(i in 1:length(Icand2fast)){ # no true, no need to fit just pred
  if(any(chosen[,2]==i)){Icand2fast[i] <- 0}else{
    Icand2fast[i] <- IMSPE2fast(g, g[i], Iselect$fit1new, Iselect$fit2new, Iselect$fit3new, constant=TRUE)$IMSPE
  }
}
for(i in 1:length(Icand3fast)){ # no true, no need to fit just pred
  if(any(chosen[,2]==i)){Icand3fast[i] <- 0}else{
    Icand3fast[i] <- IMSPE3fast(g, g[i], Iselect$fit1new, Iselect$fit2new, Iselect$fit3new, constant=TRUE)$IMSPE
  }
}
if(any(Icand1fast==0)){Icand1fast[which(Icand1fast==0)] <-  max(Icand1fast)}
if(any(Icand2fast==0)){Icand2fast[which(Icand2fast==0)] <-  max(Icand2fast)}
if(any(Icand3fast==0)){Icand3fast[which(Icand3fast==0)] <-  max(Icand3fast)}

plot(g, Icand1fast, type="l", lwd=2, col=3, ylim=range(Icand1fast))
plot(g, Icand2fast, type="l", lwd=2, col=3, ylim=range(Icand2fast))
plot(g, Icand3fast, type="l", lwd=2, col=3, ylim=range(Icand3fast))

which.min(Icand1fast)
which.min(Icand2fast)
which.min(Icand3fast)

### Fast update; Equation 6.6. in Surrogates ###
### ALC; How much can be improved. Equation 6.6. in Surrogates ###
alcfast <- c(Icurrent - Icand1fast[which.min(Icand1fast)], Icurrent - Icand2fast[which.min(Icand2fast)], Icurrent - Icand3fast[which.min(Icand3fast)] )
alcfast

### cost; 1, 2, 3 ###
which.max(alcfast/c(2,(2+8),(2+8+32)))
alcfast/c(2,(2+8),(2+8+32))


chosen <- rbind(chosen, c(which.max(alcfast/c(2,(2+8),(2+8+32))), which.min(cbind(Icand1fast, Icand2fast, Icand3fast)[,which.max(alcfast/c(2,(2+8),(2+8+32)))])))


### Plotting the chosen point ###
plot(x, predy, type="l", lwd=2, col=3, 
     ylim=range(c(predy+1.96*sqrt(predsig2*length(y3)/(length(y3)-2)), predy-1.96*sqrt(predsig2*length(y3)/(length(y3)-2))
     ))) 
lines(x, predy+1.96*sqrt(predsig2*length(y3)/(length(y3)-2)), col=3, lty=2)
lines(x, predy-1.96*sqrt(predsig2*length(y3)/(length(y3)-2)), col=3, lty=2)

curve(fl(x,l=5),add=TRUE, col=1,lwd=2,lty=2) # high fidelity(TRUE); Black

points(X1, y1, pch="1", col="red")
points(X2, y2, pch="2", col="red")
points(X3, y3, pch="3", col="red")

text(g[which.min(Icand1fast)], predy[which.min(Icand1fast)], expression("1*"), col="red")
text(g[which.min(Icand2fast)], predy[which.min(Icand2fast)], expression("2*"), col="red")
text(g[which.min(Icand3fast)], predy[which.min(Icand3fast)], expression("3*"), col="red")



#################
### Add point after 1 at 1###
#################
Iselect <- IMSPE1select(g, g[which.min(Icand1fast)], Iselect$fit1new, Iselect$fit2new, Iselect$fit3new, constant=TRUE)

### closed ###
predy <- closed2(x, Iselect$fit1new, Iselect$fit2new, Iselect$fit3new, constant=TRUE)$mu
predsig2 <- closed2(x, Iselect$fit1new, Iselect$fit2new, Iselect$fit3new, constant=TRUE)$sig2

### RMSE ###  0.07955958 -> 0.07413344 -> 0.07420309 -> 0.07138561 -> 0.07247886
sqrt(mean((predy-fl(x, l=Inf))^2)) # closed form

#############
### IMSPE ###
#############
g <- seq(0,1,0.01) # more than 1-dim, use expand.grid()
Icurrent <- mean(closed2(g, Iselect$fit1new, Iselect$fit2new, Iselect$fit3new, constant=TRUE)$sig2) 
Icurrent 
# current IMSPE, 0.002260122 -> 0.00084703 -> 0.0004360591 -> 0.0001520631 -> 0.0001400953

### Add 1 points to the low-fidelity data ###
Icand1fast <- c(rep(0, length(g))) # IMSPE candidates
Icand2fast <- c(rep(0, length(g))) # IMSPE candidates
Icand3fast <- c(rep(0, length(g))) # IMSPE candidates

for(i in 1:length(Icand1fast)){ # no true, no need to fit just pred
  if(any(chosen[,2]==i)){Icand1fast[i] <- 0}else{
    Icand1fast[i] <- IMSPE1fast(g, g[i], Iselect$fit1new, Iselect$fit2new, Iselect$fit3new, constant=TRUE)$IMSPE
  }
}
for(i in 1:length(Icand2fast)){ # no true, no need to fit just pred
  if(any(chosen[,2]==i)){Icand2fast[i] <- 0}else{
    Icand2fast[i] <- IMSPE2fast(g, g[i], Iselect$fit1new, Iselect$fit2new, Iselect$fit3new, constant=TRUE)$IMSPE
  }
}
for(i in 1:length(Icand3fast)){ # no true, no need to fit just pred
  if(any(chosen[,2]==i)){Icand3fast[i] <- 0}else{
    Icand3fast[i] <- IMSPE3fast(g, g[i], Iselect$fit1new, Iselect$fit2new, Iselect$fit3new, constant=TRUE)$IMSPE
  }
}
if(any(Icand1fast==0)){Icand1fast[which(Icand1fast==0)] <-  max(Icand1fast)}
if(any(Icand2fast==0)){Icand2fast[which(Icand2fast==0)] <-  max(Icand2fast)}
if(any(Icand3fast==0)){Icand3fast[which(Icand3fast==0)] <-  max(Icand3fast)}

plot(g, Icand1fast, type="l", lwd=2, col=3, ylim=range(Icand1fast))
plot(g, Icand2fast, type="l", lwd=2, col=3, ylim=range(Icand2fast))
plot(g, Icand3fast, type="l", lwd=2, col=3, ylim=range(Icand3fast))

which.min(Icand1fast)
which.min(Icand2fast)
which.min(Icand3fast)

### Fast update; Equation 6.6. in Surrogates ###
### ALC; How much can be improved. Equation 6.6. in Surrogates ###
alcfast <- c(Icurrent - Icand1fast[which.min(Icand1fast)], Icurrent - Icand2fast[which.min(Icand2fast)], Icurrent - Icand3fast[which.min(Icand3fast)] )
alcfast

### cost; 1, 2, 3 ###
which.max(alcfast/c(2,(2+8),(2+8+32)))
alcfast/c(2,(2+8),(2+8+32))


chosen <- rbind(chosen, c(which.max(alcfast/c(2,(2+8),(2+8+32))), which.min(cbind(Icand1fast, Icand2fast, Icand3fast)[,which.max(alcfast/c(2,(2+8),(2+8+32)))])))


### Plotting the chosen point ###
plot(x, predy, type="l", lwd=2, col=3, 
     ylim=range(c(predy+1.96*sqrt(predsig2*length(y3)/(length(y3)-2)), predy-1.96*sqrt(predsig2*length(y3)/(length(y3)-2))
     ))) 
lines(x, predy+1.96*sqrt(predsig2*length(y3)/(length(y3)-2)), col=3, lty=2)
lines(x, predy-1.96*sqrt(predsig2*length(y3)/(length(y3)-2)), col=3, lty=2)

curve(fl(x,l=5),add=TRUE, col=1,lwd=2,lty=2) # high fidelity(TRUE); Black

points(X1, y1, pch="1", col="red")
points(X2, y2, pch="2", col="red")
points(X3, y3, pch="3", col="red")

text(g[which.min(Icand1fast)], predy[which.min(Icand1fast)], expression("1*"), col="red")
text(g[which.min(Icand2fast)], predy[which.min(Icand2fast)], expression("2*"), col="red")
text(g[which.min(Icand3fast)], predy[which.min(Icand3fast)], expression("3*"), col="red")




#################
### Add point at level 2 ###
#################
Iselect <- IMSPE2select(g, g[which.min(Icand2fast)], Iselect$fit1new, Iselect$fit2new, Iselect$fit3new, constant=TRUE)

### closed ###
predy <- closed2(x, Iselect$fit1new, Iselect$fit2new, Iselect$fit3new, constant=TRUE)$mu
predsig2 <- closed2(x, Iselect$fit1new, Iselect$fit2new, Iselect$fit3new, constant=TRUE)$sig2

### RMSE ###  0.07955958 -> 0.07413344 -> 0.07420309 -> 0.07138561 -> 0.07247886
sqrt(mean((predy-fl(x, l=Inf))^2)) # closed form

#############
### IMSPE ###
#############
g <- seq(0,1,0.01) # more than 1-dim, use expand.grid()
Icurrent <- mean(closed2(g, Iselect$fit1new, Iselect$fit2new, Iselect$fit3new, constant=TRUE)$sig2) 
Icurrent
# current IMSPE, 0.002260122 -> 0.00084703 -> 0.0004360591 -> 0.0001520631 -> 0.0001400953

### Add 1 points to the low-fidelity data ###
Icand1fast <- c(rep(0, length(g))) # IMSPE candidates
Icand2fast <- c(rep(0, length(g))) # IMSPE candidates
Icand3fast <- c(rep(0, length(g))) # IMSPE candidates

for(i in 1:length(Icand1fast)){ # no true, no need to fit just pred
  if(any(chosen[,2]==i)){Icand1fast[i] <- 0}else{
    Icand1fast[i] <- IMSPE1fast(g, g[i], Iselect$fit1new, Iselect$fit2new, Iselect$fit3new, constant=TRUE)$IMSPE
  }
}
for(i in 1:length(Icand2fast)){ # no true, no need to fit just pred
  if(any(chosen[,2]==i)){Icand2fast[i] <- 0}else{
    Icand2fast[i] <- IMSPE2fast(g, g[i], Iselect$fit1new, Iselect$fit2new, Iselect$fit3new, constant=TRUE)$IMSPE
  }
}
for(i in 1:length(Icand3fast)){ # no true, no need to fit just pred
  if(any(chosen[,2]==i)){Icand3fast[i] <- 0}else{
    Icand3fast[i] <- IMSPE3fast(g, g[i], Iselect$fit1new, Iselect$fit2new, Iselect$fit3new, constant=TRUE)$IMSPE
  }
}
if(any(Icand1fast==0)){Icand1fast[which(Icand1fast==0)] <-  max(Icand1fast)}
if(any(Icand2fast==0)){Icand2fast[which(Icand2fast==0)] <-  max(Icand2fast)}
if(any(Icand3fast==0)){Icand3fast[which(Icand3fast==0)] <-  max(Icand3fast)}

plot(g, Icand1fast, type="l", lwd=2, col=3, ylim=range(Icand1fast))
plot(g, Icand2fast, type="l", lwd=2, col=3, ylim=range(Icand2fast))
plot(g, Icand3fast, type="l", lwd=2, col=3, ylim=range(Icand3fast))

which.min(Icand1fast)
which.min(Icand2fast)
which.min(Icand3fast)

### Fast update; Equation 6.6. in Surrogates ###
### ALC; How much can be improved. Equation 6.6. in Surrogates ###
alcfast <- c(Icurrent - Icand1fast[which.min(Icand1fast)], Icurrent - Icand2fast[which.min(Icand2fast)], Icurrent - Icand3fast[which.min(Icand3fast)] )
alcfast

### cost; 1, 2, 3 ###
which.max(alcfast/c(2,(2+8),(2+8+32)))
alcfast/c(2,(2+8),(2+8+32))


chosen <- rbind(chosen, c(which.max(alcfast/c(2,(2+8),(2+8+32))), which.min(cbind(Icand1fast, Icand2fast, Icand3fast)[,which.max(alcfast/c(2,(2+8),(2+8+32)))])))


### Plotting the chosen point ###
plot(x, predy, type="l", lwd=2, col=3, 
     ylim=range(c(predy+1.96*sqrt(predsig2*length(y3)/(length(y3)-2)), predy-1.96*sqrt(predsig2*length(y3)/(length(y3)-2))
     ))) 
lines(x, predy+1.96*sqrt(predsig2*length(y3)/(length(y3)-2)), col=3, lty=2)
lines(x, predy-1.96*sqrt(predsig2*length(y3)/(length(y3)-2)), col=3, lty=2)

curve(fl(x,l=5),add=TRUE, col=1,lwd=2,lty=2) # high fidelity(TRUE); Black

points(X1, y1, pch="1", col="red")
points(X2, y2, pch="2", col="red")
points(X3, y3, pch="3", col="red")

text(g[which.min(Icand1fast)], predy[which.min(Icand1fast)], expression("1*"), col="red")
text(g[which.min(Icand2fast)], predy[which.min(Icand2fast)], expression("2*"), col="red")
text(g[which.min(Icand3fast)], predy[which.min(Icand3fast)], expression("3*"), col="red")





#################
### Add point at level 3 ###
#################
Iselect <- IMSPE3select(g, g[which.min(Icand3fast)], Iselect$fit1new, Iselect$fit2new, Iselect$fit3new, constant=TRUE)

### closed ###
predy <- closed2(x, Iselect$fit1new, Iselect$fit2new, Iselect$fit3new, constant=TRUE)$mu
predsig2 <- closed2(x, Iselect$fit1new, Iselect$fit2new, Iselect$fit3new, constant=TRUE)$sig2

### RMSE ###  0.07955958 -> 0.07413344 -> 0.07420309 -> 0.07138561 -> 0.07247886
sqrt(mean((predy-fl(x, l=Inf))^2)) # closed form

#############
### IMSPE ###
#############
g <- seq(0,1,0.01) # more than 1-dim, use expand.grid()
Icurrent <- mean(closed2(g, Iselect$fit1new, Iselect$fit2new, Iselect$fit3new, constant=TRUE)$sig2) 
Icurrent
# current IMSPE, 0.002260122 -> 0.00084703 -> 0.0004360591 -> 0.0001520631 -> 0.0001400953

### Add 1 points to the low-fidelity data ###
Icand1fast <- c(rep(0, length(g))) # IMSPE candidates
Icand2fast <- c(rep(0, length(g))) # IMSPE candidates
Icand3fast <- c(rep(0, length(g))) # IMSPE candidates

for(i in 1:length(Icand1fast)){ # no true, no need to fit just pred
  if(any(chosen[,2]==i)){Icand1fast[i] <- 0}else{
    Icand1fast[i] <- IMSPE1fast(g, g[i], Iselect$fit1new, Iselect$fit2new, Iselect$fit3new, constant=TRUE)$IMSPE
  }
}
for(i in 1:length(Icand2fast)){ # no true, no need to fit just pred
  if(any(chosen[,2]==i)){Icand2fast[i] <- 0}else{
    Icand2fast[i] <- IMSPE2fast(g, g[i], Iselect$fit1new, Iselect$fit2new, Iselect$fit3new, constant=TRUE)$IMSPE
  }
}
for(i in 1:length(Icand3fast)){ # no true, no need to fit just pred
  if(any(chosen[,2]==i)){Icand3fast[i] <- 0}else{
    Icand3fast[i] <- IMSPE3fast(g, g[i], Iselect$fit1new, Iselect$fit2new, Iselect$fit3new, constant=TRUE)$IMSPE
  }
}
if(any(Icand1fast==0)){Icand1fast[which(Icand1fast==0)] <-  max(Icand1fast)}
if(any(Icand2fast==0)){Icand2fast[which(Icand2fast==0)] <-  max(Icand2fast)}
if(any(Icand3fast==0)){Icand3fast[which(Icand3fast==0)] <-  max(Icand3fast)}

plot(g, Icand1fast, type="l", lwd=2, col=3, ylim=range(Icand1fast))
plot(g, Icand2fast, type="l", lwd=2, col=3, ylim=range(Icand2fast))
plot(g, Icand3fast, type="l", lwd=2, col=3, ylim=range(Icand3fast))

which.min(Icand1fast)
which.min(Icand2fast)
which.min(Icand3fast)

### Fast update; Equation 6.6. in Surrogates ###
### ALC; How much can be improved. Equation 6.6. in Surrogates ###
alcfast <- c(Icurrent - Icand1fast[which.min(Icand1fast)], Icurrent - Icand2fast[which.min(Icand2fast)], Icurrent - Icand3fast[which.min(Icand3fast)] )
alcfast

### cost; 1, 2, 3 ###
which.max(alcfast/c(2,(2+8),(2+8+32)))
alcfast/c(2,(2+8),(2+8+32))


chosen <- rbind(chosen, c(which.max(alcfast/c(2,(2+8),(2+8+32))), which.min(cbind(Icand1fast, Icand2fast, Icand3fast)[,which.max(alcfast/c(2,(2+8),(2+8+32)))])))


### Plotting the chosen point ###
plot(x, predy, type="l", lwd=2, col=3, 
     ylim=range(c(predy+1.96*sqrt(predsig2*length(y3)/(length(y3)-2)), predy-1.96*sqrt(predsig2*length(y3)/(length(y3)-2))
     ))) 
lines(x, predy+1.96*sqrt(predsig2*length(y3)/(length(y3)-2)), col=3, lty=2)
lines(x, predy-1.96*sqrt(predsig2*length(y3)/(length(y3)-2)), col=3, lty=2)

curve(fl(x,l=5),add=TRUE, col=1,lwd=2,lty=2) # high fidelity(TRUE); Black

points(X1, y1, pch="1", col="red")
points(X2, y2, pch="2", col="red")
points(X3, y3, pch="3", col="red")

text(g[which.min(Icand1fast)], predy[which.min(Icand1fast)], expression("1*"), col="red")
text(g[which.min(Icand2fast)], predy[which.min(Icand2fast)], expression("2*"), col="red")
text(g[which.min(Icand3fast)], predy[which.min(Icand3fast)], expression("3*"), col="red")



