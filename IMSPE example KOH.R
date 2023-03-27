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
set.seed(3)
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


### KOH method ###
fit.KOH3 <- fit.KOH(X1, X2, X3, y1, y2, y3)
pred.KOH3 <- pred.KOH(fit.KOH3, x)
mx2 <- pred.KOH3$mu
koh.var2 <- pred.KOH3$sig2


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
sqrt(mean((mx2-fl(x, l=Inf))^2)) # KOH
sqrt(mean((pred3new$mu-fl(x, l=Inf))^2)) # not closed form
sqrt(mean((pred3$mu-fl(x, l=Inf))^2)) # single fidelity


#############
### IMSPE ###
#############
g <- seq(0,1,0.01) # more than 1-dim, use expand.grid()
Icurrent <- mean(koh.var2) # current IMSPE
Icurrent

### Add 1 points to the low-fidelity data ###
IcandKOH1 <- c(rep(0, length(g))) # IMSPE candidates
IcandKOH2 <- c(rep(0, length(g))) # IMSPE candidates
IcandKOH3 <- c(rep(0, length(g))) # IMSPE candidates

for(i in 1:length(IcandKOH1)){ # no true, no need to fit just pred
  IcandKOH1[i] <- IMSPEKOH1(g, g[i], fit.KOH3)
}
for(i in 1:length(IcandKOH2)){ # no true, no need to fit just pred
  IcandKOH2[i] <- IMSPEKOH2(g, g[i], fit.KOH3)
}
for(i in 1:length(IcandKOH3)){ # no true, no need to fit just pred
  IcandKOH3[i] <- IMSPEKOH3(g, g[i], fit.KOH3)
}

plot(g, IcandKOH1, type="l", lwd=2, col=3, ylim=range(IcandKOH1))
plot(g, IcandKOH2, type="l", lwd=2, col=3, ylim=range(IcandKOH2))
plot(g, IcandKOH3, type="l", lwd=2, col=3, ylim=range(IcandKOH3))

which.min(IcandKOH1)
which.min(IcandKOH2)
which.min(IcandKOH3)

### Fast update; Equation 6.6. in Surrogates ###
### ALC; How much can be improved. Equation 6.6. in Surrogates ###
alcfast <- c(Icurrent - IcandKOH1[which.min(IcandKOH1)], Icurrent - IcandKOH2[which.min(IcandKOH2)], Icurrent - IcandKOH3[which.min(IcandKOH3)] )
alcfast

### cost; 1, 2, 3 ###
which.max(alcfast/c(1,(1+2),(1+2+3)))
alcfast/c(1,(1+2),(1+2+3))

### cost; 1, 10, 100 ###
which.max(alcfast/c(1,(1+10),(1+10+100)))
alcfast/c(1,(1+10),(1+10+100))


chosen <- matrix(0, ncol=2)
chosen[1,1] <- which.max(alcfast/c(1,(1+2),(1+2+3)))
chosen[1,2] <- which.min(cbind(IcandKOH1, IcandKOH2, IcandKOH3)[,chosen[1,1]])


### Plotting the chosen point ###
plot(x, mx2, type="l", lwd=2, col=3, 
     ylim=range(c(mx2+1.96*sqrt(koh.var2*length(y3)/(length(y3)-2)), mx2-1.96*sqrt(koh.var2*length(y3)/(length(y3)-2))
     ))) 
lines(x, mx2+1.96*sqrt(koh.var2*length(y3)/(length(y3)-2)), col=3, lty=2)
lines(x, mx2-1.96*sqrt(koh.var2*length(y3)/(length(y3)-2)), col=3, lty=2)

curve(fl(x,l=5),add=TRUE, col=1,lwd=2,lty=2) # high fidelity(TRUE); Black

points(X1, y1, pch="1", col="red")
points(X2, y2, pch="2", col="red")
points(X3, y3, pch="3", col="red")

text(g[which.min(IcandKOH1)], mx2[which.min(IcandKOH1)], expression("1*"), col="red")
text(g[which.min(IcandKOH2)], mx2[which.min(IcandKOH2)], expression("2*"), col="red")
text(g[which.min(IcandKOH3)], mx2[which.min(IcandKOH3)], expression("3*"), col="red")



#################
### Add point ###
#################
Iselect <- IMSPEKOH1select(g, g[which.min(IcandKOH1)], fit.KOH3)

### closed ###
mx2 <- pred.KOH(Iselect$fit, x)$mu
koh.var2 <- pred.KOH(Iselect$fit, x)$sig2

### RMSE ###  0.07955958 -> 0.07413344 -> 0.07678651 -> 
sqrt(mean((mx2-fl(x, l=Inf))^2)) # closed form

#############
### IMSPE ###
#############
g <- seq(0,1,0.01) # more than 1-dim, use expand.grid()
Icurrent <- mean(koh.var2) # current IMSPE, 0.002260122 -> 0.00084703 -> 0.001532291 -> 
Icurrent
### Add 1 points to the low-fidelity data ###
IcandKOH1 <- c(rep(0, length(g))) # IMSPE candidates
IcandKOH2 <- c(rep(0, length(g))) # IMSPE candidates
IcandKOH3 <- c(rep(0, length(g))) # IMSPE candidates

for(i in 1:length(IcandKOH1)){ # no true, no need to fit just pred
  if(any(chosen[,2]==i)){IcandKOH1[i] <- 0}else{
    IcandKOH1[i] <- IMSPEKOH1(g, g[i], Iselect$fit)
  }
}
for(i in 1:length(IcandKOH2)){ # no true, no need to fit just pred
  if(any(chosen[,2]==i)){IcandKOH2[i] <- 0}else{
    IcandKOH2[i] <- IMSPEKOH2(g, g[i], Iselect$fit)
  }
}
for(i in 1:length(IcandKOH3)){ # no true, no need to fit just pred
  if(any(chosen[,2]==i)){IcandKOH3[i] <- 0}else{
    IcandKOH3[i] <- IMSPEKOH3(g, g[i], Iselect$fit)
  }
}
if(any(IcandKOH1==0)){IcandKOH1[which(IcandKOH1==0)] <-  max(IcandKOH1)}
if(any(IcandKOH2==0)){IcandKOH2[which(IcandKOH2==0)] <-  max(IcandKOH2)}
if(any(IcandKOH3==0)){IcandKOH3[which(IcandKOH3==0)] <-  max(IcandKOH3)}

plot(g, IcandKOH1, type="l", lwd=2, col=3, ylim=range(IcandKOH1))
plot(g, IcandKOH2, type="l", lwd=2, col=3, ylim=range(IcandKOH2))
plot(g, IcandKOH3, type="l", lwd=2, col=3, ylim=range(IcandKOH3))

which.min(IcandKOH1)
which.min(IcandKOH2)
which.min(IcandKOH3)

### Fast update; Equation 6.6. in Surrogates ###
### ALC; How much can be improved. Equation 6.6. in Surrogates ###
alcfast <- c(Icurrent - IcandKOH1[which.min(IcandKOH1)], Icurrent - IcandKOH2[which.min(IcandKOH2)], Icurrent - IcandKOH3[which.min(IcandKOH3)] )
alcfast

### cost; 1, 2, 3 ###
which.max(alcfast/c(1,(1+2),(1+2+3)))
alcfast/c(1,(1+2),(1+2+3))

### cost; 1, 10, 100 ###
which.max(alcfast/c(1,(1+10),(1+10+100)))
alcfast/c(1,(1+10),(1+10+100))


chosen <- rbind(chosen, c(which.max(alcfast/c(1,(1+2),(1+2+3))), which.min(cbind(IcandKOH1, IcandKOH2, IcandKOH3)[,which.max(alcfast/c(1,(1+2),(1+2+3)))])))


### Plotting the chosen point ###
plot(x, mx2, type="l", lwd=2, col=3, 
     ylim=range(c(mx2+1.96*sqrt(koh.var2*length(y3)/(length(y3)-2)), mx2-1.96*sqrt(koh.var2*length(y3)/(length(y3)-2))
     ))) 
lines(x, mx2+1.96*sqrt(koh.var2*length(y3)/(length(y3)-2)), col=3, lty=2)
lines(x, mx2-1.96*sqrt(koh.var2*length(y3)/(length(y3)-2)), col=3, lty=2)

curve(fl(x,l=5),add=TRUE, col=1,lwd=2,lty=2) # high fidelity(TRUE); Black

points(X1, y1, pch="1", col="red")
points(X2, y2, pch="2", col="red")
points(X3, y3, pch="3", col="red")

text(g[which.min(IcandKOH1)], mx2[which.min(IcandKOH1)], expression("1*"), col="red")
text(g[which.min(IcandKOH2)], mx2[which.min(IcandKOH2)], expression("2*"), col="red")
text(g[which.min(IcandKOH3)], mx2[which.min(IcandKOH3)], expression("3*"), col="red")



#################
### Add point after 1 at 1###
#################
Iselect <- IMSPEKOH1select(g, g[which.min(IcandKOH1)], Iselect$fit)

### closed ###
mx2 <- pred.KOH(Iselect$fit, x)$mu
koh.var2 <- pred.KOH(Iselect$fit, x)$sig2

### RMSE ###  0.07955958 -> 0.07413344 -> 0.07678651 -> 
sqrt(mean((mx2-fl(x, l=Inf))^2)) # closed form

#############
### IMSPE ###
#############
g <- seq(0,1,0.01) # more than 1-dim, use expand.grid()
Icurrent <- mean(koh.var2) 
Icurrent
# current IMSPE, 0.002260122 -> 0.00084703 -> 0.0004360591 -> 0.0001520631 -> 0.0001400953

### Add 1 points to the low-fidelity data ###
IcandKOH1 <- c(rep(0, length(g))) # IMSPE candidates
IcandKOH2 <- c(rep(0, length(g))) # IMSPE candidates
IcandKOH3 <- c(rep(0, length(g))) # IMSPE candidates

for(i in 1:length(IcandKOH1)){ # no true, no need to fit just pred
  if(any(chosen[,2]==i)){IcandKOH1[i] <- 0}else{
    IcandKOH1[i] <- IMSPEKOH1(g, g[i], Iselect$fit)
  }
}
for(i in 1:length(IcandKOH2)){ # no true, no need to fit just pred
  if(any(chosen[,2]==i)){IcandKOH2[i] <- 0}else{
    IcandKOH2[i] <- IMSPEKOH2(g, g[i], Iselect$fit)
  }
}
for(i in 1:length(IcandKOH3)){ # no true, no need to fit just pred
  if(any(chosen[,2]==i)){IcandKOH3[i] <- 0}else{
    IcandKOH3[i] <- IMSPEKOH3(g, g[i], Iselect$fit)
  }
}
if(any(IcandKOH1==0)){IcandKOH1[which(IcandKOH1==0)] <-  max(IcandKOH1)}
if(any(IcandKOH2==0)){IcandKOH2[which(IcandKOH2==0)] <-  max(IcandKOH2)}
if(any(IcandKOH3==0)){IcandKOH3[which(IcandKOH3==0)] <-  max(IcandKOH3)}

plot(g, IcandKOH1, type="l", lwd=2, col=3, ylim=range(IcandKOH1))
plot(g, IcandKOH2, type="l", lwd=2, col=3, ylim=range(IcandKOH2))
plot(g, IcandKOH3, type="l", lwd=2, col=3, ylim=range(IcandKOH3))

which.min(IcandKOH1)
which.min(IcandKOH2)
which.min(IcandKOH3)

### Fast update; Equation 6.6. in Surrogates ###
### ALC; How much can be improved. Equation 6.6. in Surrogates ###
alcfast <- c(Icurrent - IcandKOH1[which.min(IcandKOH1)], Icurrent - IcandKOH2[which.min(IcandKOH2)], Icurrent - IcandKOH3[which.min(IcandKOH3)] )
alcfast

### cost; 1, 2, 3 ###
which.max(alcfast/c(1,(1+2),(1+2+3)))
alcfast/c(1,(1+2),(1+2+3))

### cost; 1, 10, 100 ###
which.max(alcfast/c(1,(1+10),(1+10+100)))
alcfast/c(1,(1+10),(1+10+100))


chosen <- rbind(chosen, c(which.max(alcfast/c(1,(1+2),(1+2+3))), which.min(cbind(IcandKOH1, IcandKOH2, IcandKOH3)[,which.max(alcfast/c(1,(1+2),(1+2+3)))])))


### Plotting the chosen point ###
plot(x, mx2, type="l", lwd=2, col=3, 
     ylim=range(c(mx2+1.96*sqrt(koh.var2*length(y3)/(length(y3)-2)), mx2-1.96*sqrt(koh.var2*length(y3)/(length(y3)-2))
     ))) 
lines(x, mx2+1.96*sqrt(koh.var2*length(y3)/(length(y3)-2)), col=3, lty=2)
lines(x, mx2-1.96*sqrt(koh.var2*length(y3)/(length(y3)-2)), col=3, lty=2)

curve(fl(x,l=5),add=TRUE, col=1,lwd=2,lty=2) # high fidelity(TRUE); Black

points(X1, y1, pch="1", col="red")
points(X2, y2, pch="2", col="red")
points(X3, y3, pch="3", col="red")

text(g[which.min(IcandKOH1)], mx2[which.min(IcandKOH1)], expression("1*"), col="red")
text(g[which.min(IcandKOH2)], mx2[which.min(IcandKOH2)], expression("2*"), col="red")
text(g[which.min(IcandKOH3)], mx2[which.min(IcandKOH3)], expression("3*"), col="red")




#################
### Add point at level 2 ###
#################
Iselect <- IMSPEKOH2select(g, g[which.min(IcandKOH2)], Iselect$fit)

### closed ###
mx2 <- pred.KOH(Iselect$fit, x)$mu
koh.var2 <- pred.KOH(Iselect$fit, x)$sig2

### RMSE ###  0.07955958 -> 0.07413344 -> 0.07678651 -> 
sqrt(mean((mx2-fl(x, l=Inf))^2)) # closed form

#############
### IMSPE ###
#############
g <- seq(0,1,0.01) # more than 1-dim, use expand.grid()
Icurrent <- mean(koh.var2) 
Icurrent
# current IMSPE, 0.002260122 -> 0.00084703 -> 0.0004360591 -> 0.0001520631 -> 0.0001400953

### Add 1 points to the low-fidelity data ###
IcandKOH1 <- c(rep(0, length(g))) # IMSPE candidates
IcandKOH2 <- c(rep(0, length(g))) # IMSPE candidates
IcandKOH3 <- c(rep(0, length(g))) # IMSPE candidates

for(i in 1:length(IcandKOH1)){ # no true, no need to fit just pred
  if(any(chosen[,2]==i)){IcandKOH1[i] <- 0}else{
    IcandKOH1[i] <- IMSPEKOH1(g, g[i], Iselect$fit)
  }
}
for(i in 1:length(IcandKOH2)){ # no true, no need to fit just pred
  if(any(chosen[,2]==i)){IcandKOH2[i] <- 0}else{
    IcandKOH2[i] <- IMSPEKOH2(g, g[i], Iselect$fit)
  }
}
for(i in 1:length(IcandKOH3)){ # no true, no need to fit just pred
  if(any(chosen[,2]==i)){IcandKOH3[i] <- 0}else{
    IcandKOH3[i] <- IMSPEKOH3(g, g[i], Iselect$fit)
  }
}
if(any(IcandKOH1==0)){IcandKOH1[which(IcandKOH1==0)] <-  max(IcandKOH1)}
if(any(IcandKOH2==0)){IcandKOH2[which(IcandKOH2==0)] <-  max(IcandKOH2)}
if(any(IcandKOH3==0)){IcandKOH3[which(IcandKOH3==0)] <-  max(IcandKOH3)}

plot(g, IcandKOH1, type="l", lwd=2, col=3, ylim=range(IcandKOH1))
plot(g, IcandKOH2, type="l", lwd=2, col=3, ylim=range(IcandKOH2))
plot(g, IcandKOH3, type="l", lwd=2, col=3, ylim=range(IcandKOH3))

which.min(IcandKOH1)
which.min(IcandKOH2)
which.min(IcandKOH3)

### Fast update; Equation 6.6. in Surrogates ###
### ALC; How much can be improved. Equation 6.6. in Surrogates ###
alcfast <- c(Icurrent - IcandKOH1[which.min(IcandKOH1)], Icurrent - IcandKOH2[which.min(IcandKOH2)], Icurrent - IcandKOH3[which.min(IcandKOH3)] )
alcfast

### cost; 1, 2, 3 ###
which.max(alcfast/c(1,(1+2),(1+2+3)))
alcfast/c(1,(1+2),(1+2+3))

### cost; 1, 10, 100 ###
which.max(alcfast/c(1,(1+10),(1+10+100)))
alcfast/c(1,(1+10),(1+10+100))


chosen <- rbind(chosen, c(which.max(alcfast/c(1,(1+2),(1+2+3))), which.min(cbind(IcandKOH1, IcandKOH2, IcandKOH3)[,which.max(alcfast/c(1,(1+2),(1+2+3)))])))


### Plotting the chosen point ###
plot(x, mx2, type="l", lwd=2, col=3, 
     ylim=range(c(mx2+1.96*sqrt(koh.var2*length(y3)/(length(y3)-2)), mx2-1.96*sqrt(koh.var2*length(y3)/(length(y3)-2))
     ))) 
lines(x, mx2+1.96*sqrt(koh.var2*length(y3)/(length(y3)-2)), col=3, lty=2)
lines(x, mx2-1.96*sqrt(koh.var2*length(y3)/(length(y3)-2)), col=3, lty=2)

curve(fl(x,l=5),add=TRUE, col=1,lwd=2,lty=2) # high fidelity(TRUE); Black

points(X1, y1, pch="1", col="red")
points(X2, y2, pch="2", col="red")
points(X3, y3, pch="3", col="red")

text(g[which.min(IcandKOH1)], mx2[which.min(IcandKOH1)], expression("1*"), col="red")
text(g[which.min(IcandKOH2)], mx2[which.min(IcandKOH2)], expression("2*"), col="red")
text(g[which.min(IcandKOH3)], mx2[which.min(IcandKOH3)], expression("3*"), col="red")





#################
### Add point at level 3 ###
#################
Iselect <- IMSPEKOH3select(g, g[which.min(IcandKOH3)], Iselect$fit)

### closed ###
mx2 <- pred.KOH(Iselect$fit, x)$mu
koh.var2 <- pred.KOH(Iselect$fit, x)$sig2

### RMSE ###  0.07955958 -> 0.07413344 -> 0.07678651 -> 
sqrt(mean((mx2-fl(x, l=Inf))^2)) # closed form

#############
### IMSPE ###
#############
g <- seq(0,1,0.01) # more than 1-dim, use expand.grid()
Icurrent <- mean(koh.var2) 
Icurrent
# current IMSPE, 0.002260122 -> 0.00084703 -> 0.0004360591 -> 0.0001520631 -> 0.0001400953

### Add 1 points to the low-fidelity data ###
IcandKOH1 <- c(rep(0, length(g))) # IMSPE candidates
IcandKOH2 <- c(rep(0, length(g))) # IMSPE candidates
IcandKOH3 <- c(rep(0, length(g))) # IMSPE candidates

for(i in 1:length(IcandKOH1)){ # no true, no need to fit just pred
  if(any(chosen[,2]==i)){IcandKOH1[i] <- 0}else{
    IcandKOH1[i] <- IMSPEKOH1(g, g[i], Iselect$fit)
  }
}
for(i in 1:length(IcandKOH2)){ # no true, no need to fit just pred
  if(any(chosen[,2]==i)){IcandKOH2[i] <- 0}else{
    IcandKOH2[i] <- IMSPEKOH2(g, g[i], Iselect$fit)
  }
}
for(i in 1:length(IcandKOH3)){ # no true, no need to fit just pred
  if(any(chosen[,2]==i)){IcandKOH3[i] <- 0}else{
    IcandKOH3[i] <- IMSPEKOH3(g, g[i], Iselect$fit)
  }
}
if(any(IcandKOH1==0)){IcandKOH1[which(IcandKOH1==0)] <-  max(IcandKOH1)}
if(any(IcandKOH2==0)){IcandKOH2[which(IcandKOH2==0)] <-  max(IcandKOH2)}
if(any(IcandKOH3==0)){IcandKOH3[which(IcandKOH3==0)] <-  max(IcandKOH3)}

plot(g, IcandKOH1, type="l", lwd=2, col=3, ylim=range(IcandKOH1))
plot(g, IcandKOH2, type="l", lwd=2, col=3, ylim=range(IcandKOH2))
plot(g, IcandKOH3, type="l", lwd=2, col=3, ylim=range(IcandKOH3))

which.min(IcandKOH1)
which.min(IcandKOH2)
which.min(IcandKOH3)

### Fast update; Equation 6.6. in Surrogates ###
### ALC; How much can be improved. Equation 6.6. in Surrogates ###
alcfast <- c(Icurrent - IcandKOH1[which.min(IcandKOH1)], Icurrent - IcandKOH2[which.min(IcandKOH2)], Icurrent - IcandKOH3[which.min(IcandKOH3)] )
alcfast

### cost; 1, 2, 3 ###
which.max(alcfast/c(1,(1+2),(1+2+3)))
alcfast/c(1,(1+2),(1+2+3))

### cost; 1, 10, 100 ###
which.max(alcfast/c(1,(1+10),(1+10+100)))
alcfast/c(1,(1+10),(1+10+100))


chosen <- rbind(chosen, c(which.max(alcfast/c(1,(1+2),(1+2+3))), which.min(cbind(IcandKOH1, IcandKOH2, IcandKOH3)[,which.max(alcfast/c(1,(1+2),(1+2+3)))])))


### Plotting the chosen point ###
plot(x, mx2, type="l", lwd=2, col=3, 
     ylim=range(c(mx2+1.96*sqrt(koh.var2*length(y3)/(length(y3)-2)), mx2-1.96*sqrt(koh.var2*length(y3)/(length(y3)-2))
     ))) 
lines(x, mx2+1.96*sqrt(koh.var2*length(y3)/(length(y3)-2)), col=3, lty=2)
lines(x, mx2-1.96*sqrt(koh.var2*length(y3)/(length(y3)-2)), col=3, lty=2)

curve(fl(x,l=5),add=TRUE, col=1,lwd=2,lty=2) # high fidelity(TRUE); Black

points(X1, y1, pch="1", col="red")
points(X2, y2, pch="2", col="red")
points(X3, y3, pch="3", col="red")

text(g[which.min(IcandKOH1)], mx2[which.min(IcandKOH1)], expression("1*"), col="red")
text(g[which.min(IcandKOH2)], mx2[which.min(IcandKOH2)], expression("2*"), col="red")
text(g[which.min(IcandKOH3)], mx2[which.min(IcandKOH3)], expression("3*"), col="red")







