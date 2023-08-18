install.packages("devtools")
library(devtools)
install_github("cran/MuFiCokriging")
library(MuFiCokriging)
library(plgp)

library(lhs)
library(laGP)
source("GP.R")
source("closed.R")
source("IMSPE1.R")
source("IMSPE2.R")
source("IMSPE3.R")
source("/Users/junoh/Downloads/StackingDesign-Reproducibility/stacking_design.R")
library(matlabr)
source("KOH.R")

fit.KOH <- function(X1, X2, Y1, Y2, g=eps){ # need to change function for another example
  
  ### Generate output
  d1 <- data.frame(X2*0.5+0.25, rep(0.05, nrow(X2))) # scale X to [-1,1]
  write.csv(d1, "/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/temp_to_matlab.txt", row.names=F)
  run_matlab_script("/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/SolveJetBlade.m", verbose = FALSE, desktop = FALSE, 
                    splash = FALSE, display = FALSE, wait = TRUE, single_thread = FALSE,
                    intern = TRUE)
  d2 <- read.table("/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/temp_to_r.txt", sep = ",")
  Y1d2 <- d2$V4
  
  ### estimating first order ###
  fit.KOHGP1 <- KOHGP(X1, Y1)
  b1 <- 1/fit.KOHGP1$theta
  sig2_1 <- fit.KOHGP1$tau2hat
  
  ### estimating second order ###
  # KOH(X2, Y2, Y1d2)
  rho1 <- KOH(X2, Y2, Y1d2)$rho
  b2 <- 1/KOH(X2, Y2, Y1d2)$theta
  sig2_2 <- KOH(X2, Y2, Y1d2)$tau2hat
  
  return(list(b=c(b1, b2), rho=rho1, tau2hat=c(sig2_1, sig2_2), g=g, X1=X1, X2=X2, Y1=Y1, Y2=Y2))
}

pred.KOH <- function(fit, x){ # need to change function for another example
  
  X1 <- fit$X1
  X2 <- fit$X2
  Y1 <- fit$Y1
  Y2 <- fit$Y2
  
  b <- fit$b
  rho <- fit$rho
  tau2hat <- fit$tau2hat
  g <- fit$g
  
  
  ### prediction of 2nd order KOH ###
  tx1 <- cbind(rho*tau2hat[1]*covar.sep(x, X1, d=1/b[1:2], g=g), 
               rho^2*tau2hat[1]*covar.sep(x, X2, d=1/b[1:2], g=g) + tau2hat[2]*covar.sep(x, X2, d=1/b[3:4], g=g))
  
  V1 <- tau2hat[1]*covar.sep(X1, d=1/b[1:2], g=g)
  V12 <- rho*tau2hat[1]*covar.sep(X1, X2, d=1/b[1:2], g=0)
  V2 <- rho^2*tau2hat[1]*covar.sep(X2, d=1/b[1:2], g=g) + tau2hat[2]*covar.sep(X2, d=1/b[3:4], g=g)
  
  V_2 <- rbind(cbind(V1, V12), cbind(t(V12), V2))
  
  mx1 <- tx1 %*% solve(V_2) %*% c(Y1, Y2)
  
  ### posterior variance ###
  koh.var1 <- pmax(0, diag(tau2hat[2]*covar.sep(as.matrix(x), d=1/b[3:4], g=g) + tau2hat[1]*rho^2*covar.sep(as.matrix(x), d=1/b[1:2], g=g) - tx1 %*% solve(V_2)%*%t(tx1)))
  
  return(list(mu=mx1, sig2=koh.var1))
}

# costmatc <- list(NA)
# rmsematc <- list(NA)

eps <- sqrt(.Machine$double.eps)

### test data ###
d <- 2     
n <- 100
set.seed(1)
X.test <- maximinLHS(n, d)
d1 <- data.frame(X.test*0.5+0.25, rep(0.025, nrow(X.test))) # scale X to [-1,1]
write.csv(d1, "/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/temp_to_matlab.txt", row.names=F)
write.csv(X.test, "/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/temp_to_X.txt", row.names=F)
run_matlab_script("/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/SolveJetBlade.m", verbose = FALSE, desktop = FALSE, 
                  splash = FALSE, display = FALSE, wait = TRUE, single_thread = FALSE,
                  intern = TRUE)
d2 <- read.table("/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/temp_to_r.txt", sep = ",")
y.test <- d2$V4 # MM=0.3 takes 10 minutes #


### training data ###
n1 <- 10; n2 <- 5

set.seed(1) 

### Generate Input ###
X1 <- maximinLHS(n1, d)
X2 <- maximinLHS(n2, d)

NestDesign <- NestedDesignBuild(design = list(X1,X2))

X1 <- NestDesign$PX
X2 <- ExtractNestDesign(NestDesign,2)

### Y1 ###
d1 <- data.frame(X1*0.5+0.25, rep(0.05, nrow(X1))) # scale X to [-1,1]
write.csv(X1, "/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/temp_to_X.txt", row.names=F)
write.csv(d1, "/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/temp_to_matlab.txt", row.names=F)
run_matlab_script("/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/SolveJetBlade.m", verbose = FALSE, desktop = FALSE, 
                  splash = FALSE, display = FALSE, wait = TRUE, single_thread = FALSE,
                  intern = TRUE)
d2 <- read.table("/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/temp_to_r.txt", sep = ",")
y1 <- d2$V4

### Y3 ###
d1 <- data.frame(X2*0.5+0.25, rep(0.025, nrow(X2))) # scale X to [-1,1]
write.csv(X2, "/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/temp_to_X.txt", row.names=F)
write.csv(d1, "/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/temp_to_matlab.txt", row.names=F)
run_matlab_script("/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/SolveJetBlade.m", verbose = FALSE, desktop = FALSE, 
                  splash = FALSE, display = FALSE, wait = TRUE, single_thread = FALSE,
                  intern = TRUE)
d2 <- read.table("/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/temp_to_r.txt", sep = ",")
y2 <- d2$V4


### closed ###
fit.closed <- closed(X1, y1, X2, y2, constant=TRUE)
predy <- predclosed(fit.closed, X.test)$mu
predsig2 <- predclosed(fit.closed, X.test)$sig2


### model fitting for f1 ###
fit.GP1 <- GP(X1, y1, constant=TRUE)

### model fitting using (x2, f1(x2)) ###
w1.x2 <- pred.GP(fit.GP1, X2)$mu # can interpolate; nested
X2new <- cbind(X2, w1.x2) # combine (X2, f1(x2))
fit.GP2new <- GP(X2new, y2, constant=TRUE) # model fitting for f_M(X2, f1(x2))

### direct fitting ###
w1.x <- c(rep(NA, nrow(X.test)))
for(j in 1:nrow(X.test)){
  w1.x[j] <- mean(rnorm(10000, mean=pred.GP(fit.GP1, X.test)$mu[j], sd=sqrt(pred.GP(fit.GP1, X.test)$sig2[j])))
}

xxnew <- cbind(X.test, w1.x)
pred2new <- pred.GP(fit.GP2new, xxnew) # not closed form

### single fidelity ###
fit.GP2 <- GP(X2, y2, constant=TRUE)
pred2 <- pred.GP(fit.GP2, X.test) 

### KOH ###
fit.KOH2 <- fit.KOH(X1, X2, y1, y2)
pred.KOH2 <- pred.KOH(fit.KOH2, X.test)
mx2 <- pred.KOH2$mu
koh.var2 <- pred.KOH2$sig2



### RMSE ###
sqrt(sum((predy-y.test)^2))/(sqrt(sum((y.test)^2))) # closed form
sqrt(sum((pred2new$mu-y.test)^2))/(sqrt(sum((y.test)^2))) # not closed form
sqrt(sum((pred2$mu-y.test)^2))/(sqrt(sum((y.test)^2))) # single fidelity
sqrt(sum((mx2-y.test)^2))/(sqrt(sum((y.test)^2))) # KOH


### IMSPE ###
Icurrent <- mean(koh.var2) # current IMSPE
Icurrent

### Add 1 points and calculate IMSPE ###
IcandKOH1 <- c(rep(0, nrow(X.test))) # IMSPE candidates
IcandKOH2 <- c(rep(0, nrow(X.test))) # IMSPE candidates

for(i in 1:length(IcandKOH1)){ # no true, no need to fit just pred
  IcandKOH1[i] <- IMSPEKOH(X.test, X.test[i,], fit.KOH2, level=1)
}
for(i in 1:length(IcandKOH2)){ # no true, no need to fit just pred
  IcandKOH2[i] <- IMSPEKOH(X.test, X.test[i,], fit.KOH2, level=2)
}

which.min(IcandKOH1)
which.min(IcandKOH2)

### Fast update; Equation 6.6. in Surrogates ###
### ALC; How much can be improved. Equation 6.6. in Surrogates ###
alcfast <- c(Icurrent - IcandKOH1[which.min(IcandKOH1)], Icurrent - IcandKOH2[which.min(IcandKOH2)])
alcfast

### cost; 1, 2, 3 ###
which.max(alcfast/c(0.62, (0.62+1.03)))
alcfast/c(0.62, (0.62+1.03))


chosen <- matrix(0, ncol=2)
chosen[1,1] <- which.max(alcfast/c(0.62, (0.62+1.03)))
chosen[1,2] <- which.min(cbind(IcandKOH1, IcandKOH2)[,chosen[1,1]])



nonlinear.cost <- 0
nonlinear.error <- sqrt(sum((mx2-y.test)^2))/(sqrt(sum((y.test)^2)))

Iselect <- IMSPEKOHselect(X.test, X.test[chosen[nrow(chosen),2],], fit.KOH2, level=chosen[nrow(chosen),1]) 



#################
### Add point ###
#################
while(nonlinear.cost[length(nonlinear.cost)] < 30){ # if total cost is less than the budget
  
  ### KOH ###
  pred.KOH2 <- pred.KOH(Iselect, X.test)
  mx2 <- pred.KOH2$mu
  koh.var2 <- pred.KOH2$sig2
  
  
  ### RMSE ###  
  nonlinear.error <- c(nonlinear.error, sqrt(sum((mx2-y.test)^2))/(sqrt(sum((y.test)^2)))) # closed form
  if(chosen[nrow(chosen),1] == 1){
    nonlinear.cost[length(nonlinear.cost)+1] <- nonlinear.cost[length(nonlinear.cost)]+0.62
  }else if(chosen[nrow(chosen),1] == 2){
    nonlinear.cost[length(nonlinear.cost)+1] <- nonlinear.cost[length(nonlinear.cost)]+(0.62+1.03)
  }
  
  #############
  ### IMSPE ###
  #############
  Icurrent <- mean(koh.var2)
  
  ### Add 1 points to the low-fidelity data ###
  IcandKOH1 <- c(rep(0, nrow(X.test))) # IMSPE candidates
  IcandKOH2 <- c(rep(0, nrow(X.test))) # IMSPE candidates
  
  for(i in 1:length(IcandKOH1)){ # no true, no need to fit just pred
    if(any(chosen[,2]==i)){IcandKOH1[i] <- 0}else{
      IcandKOH1[i] <- IMSPEKOH(X.test, X.test[i,], fit.KOH2, level=1)
    }
  }
  for(i in 1:length(IcandKOH2)){ # no true, no need to fit just pred
    if(any(chosen[,2]==i)){IcandKOH2[i] <- 0}else{
      IcandKOH2[i] <- IMSPEKOH(X.test, X.test[i,], fit.KOH2, level=2)
    }
  }
  
  if(any(IcandKOH1==0)){IcandKOH1[which(IcandKOH1==0)] <-  max(IcandKOH1)}
  if(any(IcandKOH2==0)){IcandKOH2[which(IcandKOH2==0)] <-  max(IcandKOH2)}
  
  which.min(IcandKOH1)
  which.min(IcandKOH2)
  
  ### Fast update; Equation 6.6. in Surrogates ###
  ### ALC; How much can be improved. Equation 6.6. in Surrogates ###
  alcfast <- c(Icurrent - IcandKOH1[which.min(IcandKOH1)], Icurrent - IcandKOH2[which.min(IcandKOH2)])
  alcfast
  
  ### cost; 1, 2, 3 ###
  which.max(alcfast/c(0.62, (0.62+1.03)))
  alcfast/c(0.62, (0.62+1.03))
  
  
  chosen <- rbind(chosen, c(which.max(alcfast/c(0.62, (0.62+1.03))), which.min(cbind(IcandKOH1, IcandKOH2)[,which.max(alcfast/c(0.62, (0.62+1.03)))])))
  Iselect <- IMSPEKOHselect(X.test, X.test[chosen[nrow(chosen),2],], Iselect, level=chosen[nrow(chosen),1])
  
  if(nonlinear.cost[length(nonlinear.cost)] >= 100){break}
  
}


### Save results ###
nonlinear.cost
nonlinear.error


> nonlinear.cost
[1]  0.00  0.62  2.27  3.92  5.57  7.22  8.87
[8] 10.52 12.17 13.82 14.44 16.09 17.74 19.39
[15] 21.04 22.69 24.34 25.99
> nonlinear.error
[1] 0.12425586 0.11967746 0.12460465
[4] 0.09712361 0.09959276 0.11412961
[7] 0.12883777 0.12518020 0.12946475
[10] 0.20889560 0.21012229 0.20950560
[13] 0.17168923 0.16963701 0.12155827
[16] 0.12873649 0.15422820 0.14776418

