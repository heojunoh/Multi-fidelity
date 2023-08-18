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
n1 <- 12; n2 <- 8; n3 <- 5

for(kk in 1:10){
  set.seed(kk)
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

### test data ###
x <- seq(0,1,0.01)


### KOH method ###
fit.KOH3 <- fit.KOH(X1, X2, X3, y1, y2, y3)
pred.KOH3 <- pred.KOH(fit.KOH3, x)
mx2 <- pred.KOH3$mu
koh.var2 <- pred.KOH3$sig2

### RMSE ###
sqrt(mean((mx2-fl(x, l=5))^2)) # KOH

#############
### IMSPE ###
#############
### IMSPE ###
g <- seq(0,1,0.01) # more than 1-dim, use expand.grid()
Icurrent <- mean(koh.var2) # current IMSPE
Icurrent

### Add 1 points and calculate IMSPE ###
IcandKOH1 <- c(rep(0, length(g))) # IMSPE candidates
IcandKOH2 <- c(rep(0, length(g))) # IMSPE candidates

for(i in 1:length(IcandKOH1)){ # no true, no need to fit just pred
  IcandKOH1[i] <- IMSPEKOH1(g, g[i], fit.KOH2, level=1)
}
for(i in 1:length(IcandKOH2)){ # no true, no need to fit just pred
  IcandKOH2[i] <- IMSPEKOH1(g, g[i], fit.KOH2, level=2)
}

which.min(IcandKOH1)
which.min(IcandKOH2)

### Fast update; Equation 6.6. in Surrogates ###
### ALC; How much can be improved. Equation 6.6. in Surrogates ###
alcfast <- c(Icurrent - IcandKOH1[which.min(IcandKOH1)], Icurrent - IcandKOH2[which.min(IcandKOH2)], Icurrent - IcandKOH3[which.min(IcandKOH3)] )
alcfast

### cost; 2, 4, 8 ###
which.max(alcfast/c(2,4,8))
alcfast/c(2,4,8)


chosen <- matrix(0, ncol=2)
chosen[1,1] <- which.max(alcfast/c(2,4,8))
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
Iselect <- IMSPEKOHselect(g, g[which.min(IcandKOH1)], fit.KOH3, level=1)

### closed ###
mx2 <- pred.KOH(Iselect$fit, x)$mu
koh.var2 <- pred.KOH(Iselect$fit, x)$sig2

### RMSE ###  0.07955958 -> 0.07413344 -> 0.07678651 -> 
sqrt(mean((mx2-fl(x, l=5))^2)) # closed form

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
  if(any(matrix(chosen[chosen[,1]==1,], ncol=2)[,2]==i)){IcandKOH1[i] <- 0}else{
    IcandKOH1[i] <- IMSPEKOH(g, g[i], Iselect$fit, level=1)
  }
}
for(i in 1:length(IcandKOH2)){ # no true, no need to fit just pred
  if(any(matrix(chosen[chosen[,1]==2,], ncol=2)[,2]==i)){IcandKOH2[i] <- 0}else{
    IcandKOH2[i] <- IMSPEKOH(g, g[i], Iselect$fit, level=2)
  }
}
for(i in 1:length(IcandKOH3)){ # no true, no need to fit just pred
  if(any(matrix(chosen[chosen[,1]==3,], ncol=2)[,2]==i)){IcandKOH3[i] <- 0}else{
    IcandKOH3[i] <- IMSPEKOH(g, g[i], Iselect$fit, level=3)
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

### cost; 2, 4, 8 ###
which.max(alcfast/c(2,4,8))
alcfast/c(2,4,8)


chosen <- rbind(chosen, c(which.max(alcfast/c(2,4,8)), which.min(cbind(IcandKOH1, IcandKOH2, IcandKOH3)[,which.max(alcfast/c(2,4,8))])))


### Plotting the chosen point ###
plot(x, mx2, type="l", lwd=2, col=3, 
     ylim=range(c(mx2+1.96*sqrt(koh.var2*length(y3)/(length(y3)-2)), mx2-1.96*sqrt(koh.var2*length(y3)/(length(y3)-2))
     ))) 
lines(x, mx2+1.96*sqrt(koh.var2*length(y3)/(length(y3)-2)), col=3, lty=2)
lines(x, mx2-1.96*sqrt(koh.var2*length(y3)/(length(y3)-2)), col=3, lty=2)

curve(fl(x,l=5),add=TRUE, col=1,lwd=2,lty=2) # high fidelity(TRUE); Black

points(Iselect$fit$X1, Iselect$fit$Y1, pch="1", col="red")
points(Iselect$fit$X2, Iselect$fit$Y2, pch="2", col="red")
points(Iselect$fit$X3, Iselect$fit$Y3, pch="3", col="red")

text(g[which.min(IcandKOH1)], mx2[which.min(IcandKOH1)], expression("1*"), col="red")
text(g[which.min(IcandKOH2)], mx2[which.min(IcandKOH2)], expression("2*"), col="red")
text(g[which.min(IcandKOH3)], mx2[which.min(IcandKOH3)], expression("3*"), col="red")



#################
### Add point after 1 at 1###
#################
Iselect <- IMSPEKOHselect(g, g[which.min(IcandKOH1)], Iselect$fit, level=1)

### closed ###
mx2 <- pred.KOH(Iselect$fit, x)$mu
koh.var2 <- pred.KOH(Iselect$fit, x)$sig2

### RMSE ###  0.07955958 -> 0.07413344 -> 0.07678651 -> 
sqrt(mean((mx2-fl(x, l=5))^2)) # closed form

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
  if(any(matrix(chosen[chosen[,1]==1,], ncol=2)[,2]==i)){IcandKOH1[i] <- 0}else{
    IcandKOH1[i] <- IMSPEKOH(g, g[i], Iselect$fit, level=1)
  }
}
for(i in 1:length(IcandKOH2)){ # no true, no need to fit just pred
  if(any(matrix(chosen[chosen[,1]==2,], ncol=2)[,2]==i)){IcandKOH2[i] <- 0}else{
    IcandKOH2[i] <- IMSPEKOH(g, g[i], Iselect$fit, level=2)
  }
}
for(i in 1:length(IcandKOH3)){ # no true, no need to fit just pred
  if(any(matrix(chosen[chosen[,1]==3,], ncol=2)[,2]==i)){IcandKOH3[i] <- 0}else{
    IcandKOH3[i] <- IMSPEKOH(g, g[i], Iselect$fit, level=3)
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

### cost; 2, 4, 8 ###
which.max(alcfast/c(2,4,8))
alcfast/c(2,4,8)


chosen <- rbind(chosen, c(which.max(alcfast/c(2,4,8)), which.min(cbind(IcandKOH1, IcandKOH2, IcandKOH3)[,which.max(alcfast/c(2,4,8))])))


### Plotting the chosen point ###
plot(x, mx2, type="l", lwd=2, col=3, 
     ylim=range(c(mx2+1.96*sqrt(koh.var2*length(y3)/(length(y3)-2)), mx2-1.96*sqrt(koh.var2*length(y3)/(length(y3)-2))
     ))) 
lines(x, mx2+1.96*sqrt(koh.var2*length(y3)/(length(y3)-2)), col=3, lty=2)
lines(x, mx2-1.96*sqrt(koh.var2*length(y3)/(length(y3)-2)), col=3, lty=2)

curve(fl(x,l=5),add=TRUE, col=1,lwd=2,lty=2) # high fidelity(TRUE); Black

points(Iselect$fit$X1, Iselect$fit$Y1, pch="1", col="red")
points(Iselect$fit$X2, Iselect$fit$Y2, pch="2", col="red")
points(Iselect$fit$X3, Iselect$fit$Y3, pch="3", col="red")

text(g[which.min(IcandKOH1)], mx2[which.min(IcandKOH1)], expression("1*"), col="red")
text(g[which.min(IcandKOH2)], mx2[which.min(IcandKOH2)], expression("2*"), col="red")
text(g[which.min(IcandKOH3)], mx2[which.min(IcandKOH3)], expression("3*"), col="red")




#################
### Add point at level 2 ###
#################
Iselect <- IMSPEKOHselect(g, g[which.min(IcandKOH2)], Iselect$fit, level=2)

### closed ###
mx2 <- pred.KOH(Iselect$fit, x)$mu
koh.var2 <- pred.KOH(Iselect$fit, x)$sig2

### RMSE ###  0.07955958 -> 0.07413344 -> 0.07678651 -> 
sqrt(mean((mx2-fl(x, l=5))^2)) # closed form

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
  if(any(matrix(chosen[chosen[,1]==1,], ncol=2)[,2]==i)){IcandKOH1[i] <- 0}else{
    IcandKOH1[i] <- IMSPEKOH(g, g[i], Iselect$fit, level=1)
  }
}
for(i in 1:length(IcandKOH2)){ # no true, no need to fit just pred
  if(any(matrix(chosen[chosen[,1]==2,], ncol=2)[,2]==i)){IcandKOH2[i] <- 0}else{
    IcandKOH2[i] <- IMSPEKOH(g, g[i], Iselect$fit, level=2)
  }
}
for(i in 1:length(IcandKOH3)){ # no true, no need to fit just pred
  if(any(matrix(chosen[chosen[,1]==3,], ncol=2)[,2]==i)){IcandKOH3[i] <- 0}else{
    IcandKOH3[i] <- IMSPEKOH(g, g[i], Iselect$fit, level=3)
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

### cost; 2, 4, 8 ###
which.max(alcfast/c(2,4,8))
alcfast/c(2,4,8)


chosen <- rbind(chosen, c(which.max(alcfast/c(2,4,8)), which.min(cbind(IcandKOH1, IcandKOH2, IcandKOH3)[,which.max(alcfast/c(2,4,8))])))


### Plotting the chosen point ###
plot(x, mx2, type="l", lwd=2, col=3, 
     ylim=range(c(mx2+1.96*sqrt(koh.var2*length(y3)/(length(y3)-2)), mx2-1.96*sqrt(koh.var2*length(y3)/(length(y3)-2))
     ))) 
lines(x, mx2+1.96*sqrt(koh.var2*length(y3)/(length(y3)-2)), col=3, lty=2)
lines(x, mx2-1.96*sqrt(koh.var2*length(y3)/(length(y3)-2)), col=3, lty=2)

curve(fl(x,l=5),add=TRUE, col=1,lwd=2,lty=2) # high fidelity(TRUE); Black

points(Iselect$fit$X1, Iselect$fit$Y1, pch="1", col="red")
points(Iselect$fit$X2, Iselect$fit$Y2, pch="2", col="red")
points(Iselect$fit$X3, Iselect$fit$Y3, pch="3", col="red")

text(g[which.min(IcandKOH1)], mx2[which.min(IcandKOH1)], expression("1*"), col="red")
text(g[which.min(IcandKOH2)], mx2[which.min(IcandKOH2)], expression("2*"), col="red")
text(g[which.min(IcandKOH3)], mx2[which.min(IcandKOH3)], expression("3*"), col="red")





#################
### Add point at level 3 ###
#################
Iselect <- IMSPEKOHselect(g, g[which.min(IcandKOH3)], Iselect$fit, level=3)

### closed ###
mx2 <- pred.KOH(Iselect$fit, x)$mu
koh.var2 <- pred.KOH(Iselect$fit, x)$sig2

### RMSE ###  0.07955958 -> 0.07413344 -> 0.07678651 -> 
sqrt(mean((mx2-fl(x, l=5))^2)) # closed form

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
  if(any(matrix(chosen[chosen[,1]==1,], ncol=2)[,2]==i)){IcandKOH1[i] <- 0}else{
    IcandKOH1[i] <- IMSPEKOH(g, g[i], Iselect$fit, level=1)
  }
}
for(i in 1:length(IcandKOH2)){ # no true, no need to fit just pred
  if(any(matrix(chosen[chosen[,1]==2,], ncol=2)[,2]==i)){IcandKOH2[i] <- 0}else{
    IcandKOH2[i] <- IMSPEKOH(g, g[i], Iselect$fit, level=2)
  }
}
for(i in 1:length(IcandKOH3)){ # no true, no need to fit just pred
  if(any(matrix(chosen[chosen[,1]==3,], ncol=2)[,2]==i)){IcandKOH3[i] <- 0}else{
    IcandKOH3[i] <- IMSPEKOH(g, g[i], Iselect$fit, level=3)
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

### cost; 2, 4, 8 ###
which.max(alcfast/c(2,4,8))
alcfast/c(2,4,8)


chosen <- rbind(chosen, c(which.max(alcfast/c(2,4,8)), which.min(cbind(IcandKOH1, IcandKOH2, IcandKOH3)[,which.max(alcfast/c(2,4,8))])))


### Plotting the chosen point ###
plot(x, mx2, type="l", lwd=2, col=3, 
     ylim=range(c(mx2+1.96*sqrt(koh.var2*length(y3)/(length(y3)-2)), mx2-1.96*sqrt(koh.var2*length(y3)/(length(y3)-2))
     ))) 
lines(x, mx2+1.96*sqrt(koh.var2*length(y3)/(length(y3)-2)), col=3, lty=2)
lines(x, mx2-1.96*sqrt(koh.var2*length(y3)/(length(y3)-2)), col=3, lty=2)

curve(fl(x,l=5),add=TRUE, col=1,lwd=2,lty=2) # high fidelity(TRUE); Black

points(Iselect$fit$X1, Iselect$fit$Y1, pch="1", col="red")
points(Iselect$fit$X2, Iselect$fit$Y2, pch="2", col="red")
points(Iselect$fit$X3, Iselect$fit$Y3, pch="3", col="red")

text(g[which.min(IcandKOH1)], mx2[which.min(IcandKOH1)], expression("1*"), col="red")
text(g[which.min(IcandKOH2)], mx2[which.min(IcandKOH2)], expression("2*"), col="red")
text(g[which.min(IcandKOH3)], mx2[which.min(IcandKOH3)], expression("3*"), col="red")







