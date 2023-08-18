library(lhs)
library(laGP)
source("GP.R")
source("KOH.R")
source("closed.R")
source("cokm.R")

costmatco <- list(NA)
rmsematco <- list(NA)
### synthetic function ###
f1 <- function(x)
{
  sin(8*pi*x)
}

f2 <- function(x)
{ 
  (x-sqrt(2))*(sin(8*pi*x))^2
}

### training data ###
n1 <- 12; n2 <- 9


set.seed(1)
X1 <- maximinLHS(n1, 1)
X2 <- maximinLHS(n2, 1)

NestDesign <- NestedDesignBuild(design = list(X1,X2))

X1 <- NestDesign$PX
X2 <- ExtractNestDesign(NestDesign,2)

y1 <- f1(X1)
y2 <- f2(X2)

### model fitting for f1 ###
eps <- sqrt(.Machine$double.eps)
fit.GP1 <- GP(X1, y1, constant=TRUE)

### model fitting using (x2, f1(x2)) ###
w1.x2 <- pred.GP(fit.GP1, X2)$mu # can interpolate; nested
X2new <- cbind(X2, w1.x2) # combine (X2, f1(x2))
fit.GP2new <- GP(X2new, y2, constant=TRUE) # model fitting for f_M(X2, f1(x2))

### test data ###
x <- seq(0,1,0.01)


### closed ###
fit.closed <- closed(X1, y1, X2, y2, constant=TRUE)
predy <- predclosed(fit.closed, x)$mu
predsig2 <- predclosed(fit.closed, x)$sig2


### cokm ###
fit.muficokm <- MuFicokm(formula = list(~1,~1), MuFidesign = NestDesign, covtype="gauss",
                         # coef.trend = list(0,c(0,0)),
                         lower=eps, upper=0.1,
                         response = list(y1,y2), nlevel = 2)
pred.muficokm <- predict(fit.muficokm, x, "SK", cov.compute=TRUE)


### compared to single fidelity ###
fit.GP2 <- GP(X2, y2, constant=TRUE)
pred2 <- pred.GP(fit.GP2, x)

### direct fitting; not using closed form. f1(u) and f_M(u) from (u, f_M(u, f1(u))) are random variables.
w1.x <- c(rep(NA, length(x)))
for(i in 1:length(x)){
  w1.x[i] <- mean(rnorm(10000, mean=pred.GP(fit.GP1, x[i])$mu, sd=sqrt(pred.GP(fit.GP1, x[i])$sig2)))
}

xxnew <- cbind(x, w1.x)
pred2new <- pred.GP(fit.GP2new, xxnew) # not closed form


### RMSE ###
sqrt(sum((predy-f2(x))^2))/(sqrt(sum((f2(x))^2))) # closed form
sqrt(sum((pred2new$mu-f2(x))^2))/(sqrt(sum((f2(x))^2))) # not closed form
sqrt(sum((pred2$mu-f2(x))^2))/(sqrt(sum((f2(x))^2))) # single fidelity
sqrt(sum((pred.muficokm$mean-f2(x))^2))/(sqrt(sum((f2(x))^2))) # Cokm


### IMSPE ###
g <- seq(0,1,0.01) # more than 1-dim, use expand.grid()
Icurrent <- mean(pred.muficokm$sig2) # current IMSPE
Icurrent

### Add 1 points and calculate IMSPE ###
Icandcokm1 <- c(rep(0, length(g))) # IMSPE candidates
Icandcokm2 <- c(rep(0, length(g))) # IMSPE candidates

for(i in 1:length(Icandcokm1)){ # no true, no need to fit just pred
  IMSPE <- IMSPEcokm1(g, g[i], fit.muficokm)
  Icandcokm1[i] <- IMSPE$k1
  Icandcokm2[i] <- IMSPE$k2
}
# plot(g, Icandcokm1, type="l", lwd=2, col=3, ylim=range(Icandcokm1))
# plot(g, Icandcokm2, type="l", lwd=2, col=3, ylim=range(Icandcokm2))

which.min(Icandcokm1)
which.min(Icandcokm2)

### Fast update; Equation 6.6. in Surrogates ###
### ALC; How much can be improved. Equation 6.6. in Surrogates ###
alcfast <- c(Icandcokm1[which.min(Icandcokm1)], Icandcokm2[which.min(Icandcokm2)])
alcfast

### cost; 1, 2, 3 ###
which.min(alcfast*c(2,2+8))
alcfast*c(2,2+8)


chosen <- matrix(0, ncol=2)
chosen[1,1] <- which.min(alcfast*c(2,2+8))
chosen[1,2] <- which.min(cbind(Icandcokm1, Icandcokm2)[,chosen[1,1]])


### Plotting the chosen point ###
plot(x, pred.muficokm$mean, type="l", lwd=2, col=3, ylim=c(-2,1))
lines(x, pred.muficokm$mean+1.96*sqrt(pred.muficokm$sig2*length(y2)/(length(y2)-2)), col=3, lty=2)
lines(x, pred.muficokm$mean-1.96*sqrt(pred.muficokm$sig2*length(y2)/(length(y2)-2)), col=3, lty=2)

curve(f2(x),add=TRUE, col=1,lwd=2,lty=2) # high fidelity(TRUE); Black

points(X1, y1, pch="1", col="red")
points(X2, y2, pch="2", col="red")

text(g[which.min(Icandcokm1)], pred.muficokm$mean[which.min(Icandcokm1)], expression("1*"), col="red")
text(g[which.min(Icandcokm2)], pred.muficokm$mean[which.min(Icandcokm2)], expression("2*"), col="red")

nonlinear.cost <- 0
nonlinear.error <- sqrt(sum((pred.muficokm$mean-f2(x))^2))/(sqrt(sum((f2(x))^2)))

Iselect <- IMSPEcokmselect1(g, g[chosen[nrow(chosen),2]], fit.muficokm, level=chosen[nrow(chosen),1]) 


#################
### Add point ###
#################
while(nonlinear.cost[length(nonlinear.cost)] < 100){ # if total cost is less than the budget
  
  ### predictive ###
  pred.muficokm <- predict(Iselect, x, "SK", cov.compute=TRUE)
  
  ### RMSE ###  
  nonlinear.error <- c(nonlinear.error, sqrt(sum((pred.muficokm$mean-f2(x))^2))/(sqrt(sum((f2(x))^2)))) # Cokm) # closed form
  if(chosen[nrow(chosen),1] == 1){
    nonlinear.cost[length(nonlinear.cost)+1] <- nonlinear.cost[length(nonlinear.cost)]+2
  }else{
    nonlinear.cost[length(nonlinear.cost)+1] <- nonlinear.cost[length(nonlinear.cost)]+(2+8)
  }
  
  #############
  ### IMSPE ###
  #############
  g <- seq(0,1,0.01) # more than 1-dim, use expand.grid()
  Icurrent <- mean(pred.muficokm$sig2)
  
  ### Add 1 points and calculate IMSPE ###
  Icandcokm1 <- c(rep(0, length(g))) # IMSPE candidates
  Icandcokm2 <- c(rep(0, length(g))) # IMSPE candidates
  
  for(i in 1:length(Icandcokm1)){ # no true, no need to fit just pred
    if(any(chosen[,2]==i)){Icandcokm1[i] <- Icandcokm2[i] <- 0}else{
      IMSPE <- IMSPEcokm1(g, g[i], Iselect)
      Icandcokm1[i] <- IMSPE$k1
      Icandcokm2[i] <- IMSPE$k2
    }
  }
  
  if(any(Icandcokm1==0)){Icandcokm1[which(Icandcokm1==0)] <-  max(Icandcokm1)}
  if(any(Icandcokm2==0)){Icandcokm2[which(Icandcokm2==0)] <-  max(Icandcokm2)}
  
  which.min(Icandcokm1)
  which.min(Icandcokm2)

  ### Fast update; Equation 6.6. in Surrogates ###
  ### ALC; How much can be improved. Equation 6.6. in Surrogates ###
  alcfast <- c(Icandcokm1[which.min(Icandcokm1)], Icandcokm2[which.min(Icandcokm2)])
  alcfast
  
  ### cost; 1, 2, 3 ###
  which.min(alcfast*c(2,2+8))
  alcfast*c(2,2+8)
  
  
  chosen <- rbind(chosen, c(which.min(alcfast*c(2,2+8)), which.min(cbind(Icandcokm1, Icandcokm2)[,which.min(alcfast*c(2,2+8))])))
  Iselect <- IMSPEcokmselect1(g, g[chosen[nrow(chosen),2]], Iselect, level=chosen[nrow(chosen),1])
  
  if(nonlinear.cost[length(nonlinear.cost)] >= 100){break}
  
}


### Plotting the chosen point ###
plot(x, pred.muficokm$mean, type="l", lwd=2, col=3,
     ylim=c(-2,1))
lines(x, pred.muficokm$mean+1.96*sqrt(pred.muficokm$sig2*length(y2)/(length(y2)-2)), col=3, lty=2)
lines(x, pred.muficokm$mean-1.96*sqrt(pred.muficokm$sig2*length(y2)/(length(y2)-2)), col=3, lty=2)

curve(f2(x),add=TRUE, col=1,lwd=2,lty=2) # high fidelity(TRUE); Black

points(Iselect$cok[[1]]@X, Iselect$cok[[1]]@y, pch="1", col="red") 
points(Iselect$cok[[2]]@X, Iselect$cok[[2]]@y, pch="2", col="red")

text(g[which.min(Icandcokm1)], pred.muficokm$mean[which.min(Icandcokm1)], expression("1*"), col="red")
text(g[which.min(Icandcokm2)], pred.muficokm$mean[which.min(Icandcokm2)], expression("2*"), col="red")


### Save results ###
costmatco[[10]] <- nonlinear.cost
rmsematco[[10]] <- nonlinear.error
costmatco
rmsematco



