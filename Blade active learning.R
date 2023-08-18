install.packages("devtools")
library(devtools)
install_github("cran/MuFiCokriging")
library(MuFiCokriging)
library(plgp)

library(lhs)
library(laGP)
source("closed.R")
source("IMSPE1.R")
source("IMSPE2.R")
source("IMSPE3.R")
source("/Users/junoh/Downloads/StackingDesign-Reproducibility/stacking_design.R")
library(matlabr)
source("GP.R")

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

### Y2 ###
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


### RMSE ###
sqrt(sum((predy-y.test)^2))/(sqrt(sum((y.test)^2))) # closed form
sqrt(sum((pred2new$mu-y.test)^2))/(sqrt(sum((y.test)^2))) # not closed form
sqrt(sum((pred2$mu-y.test)^2))/(sqrt(sum((y.test)^2))) # single fidelity


### IMSPE ###
Icurrent <- mean(predsig2) # current IMSPE
Icurrent

### Add 1 points and calculate IMSPE ###
Icand1fast <- c(rep(0, nrow(X.test))) # IMSPE candidates
Icand2fast <- c(rep(0, nrow(X.test))) # IMSPE candidates

for(i in 1:length(Icand1fast)){ # no true, no need to fit just pred
  Icand1fast[i] <- IMSPEfast2(X.test, X.test[i,], fit.closed, level=1)$IMSPE
}
for(i in 1:length(Icand2fast)){ # no true, no need to fit just pred
  Icand2fast[i] <- IMSPEfast2(X.test, X.test[i,], fit.closed, level=2)$IMSPE
}


# which.max(predsig2)
# Icand1fast <- IMSPEfast2(X.test, X.test[which.max(predsig2),], fit.closed, level=1)$IMSPE
# Icand2fast <- IMSPEfast2(X.test, X.test[which.max(predsig2),], fit.closed, level=2)$IMSPE

### Fast update; Equation 6.6. in Surrogates ###
### ALC; How much can be improved. Equation 6.6. in Surrogates ###
alcfast <- c(Icand1fast[which.min(Icand1fast)], Icand2fast[which.min(Icand2fast)])
alcfast

### cost; 1, 2, 3 ###
which.min(alcfast*c(0.62, (0.62+1.03)))
alcfast*c(0.62, (0.62+1.03))


chosen <- matrix(0, ncol=2)
chosen[1,1] <- which.min(alcfast*c(0.62, (0.62+1.03)))
chosen[1,2] <- which.min(cbind(Icand1fast, Icand2fast)[,chosen[1,1]])



nonlinear.cost <- 0
nonlinear.error <- sqrt(sum((predy-y.test)^2))/(sqrt(sum((y.test)^2)))

Iselect <- IMSPEselect2(X.test, X.test[chosen[nrow(chosen),2],], fit.closed, level=chosen[nrow(chosen),1]) 



#################
### Add point ###
#################
while(nonlinear.cost[length(nonlinear.cost)] < 30){ # if total cost is less than the budget
  
  ### closed ###
  predy <- predclosed(Iselect$fit, X.test)$mu
  predsig2 <- predclosed(Iselect$fit, X.test)$sig2
  
  ### RMSE ###  
  nonlinear.error <- c(nonlinear.error, sqrt(sum((predy-y.test)^2))/(sqrt(sum((y.test)^2)))) # closed form
  if(chosen[nrow(chosen),1] == 1){
    nonlinear.cost[length(nonlinear.cost)+1] <- nonlinear.cost[length(nonlinear.cost)]+0.62
  }else if(chosen[nrow(chosen),1] == 2){
    nonlinear.cost[length(nonlinear.cost)+1] <- nonlinear.cost[length(nonlinear.cost)]+(0.62+1.03)
  }
  
  #############
  ### IMSPE ###
  #############
  Icurrent <- Iselect$IMSPE
  
  ### Add 1 points to the low-fidelity data ###
  Icand1fast <- c(rep(0, nrow(X.test))) # IMSPE candidates
  Icand2fast <- c(rep(0, nrow(X.test))) # IMSPE candidates
  
  for(i in 1:length(Icand1fast)){ # no true, no need to fit just pred
    if(any(chosen[,2]==i)){Icand1fast[i] <- 0}else{
      Icand1fast[i] <- IMSPEfast2(X.test, X.test[i,], fit.closed, level=1)$IMSPE
    }
  }
  for(i in 1:length(Icand2fast)){ # no true, no need to fit just pred
    if(any(chosen[,2]==i)){Icand2fast[i] <- 0}else{
      Icand2fast[i] <- IMSPEfast2(X.test, X.test[i,], fit.closed, level=2)$IMSPE
    }
  }
  
  if(any(Icand1fast==0)){Icand1fast[which(Icand1fast==0)] <-  max(Icand1fast)}
  if(any(Icand2fast==0)){Icand2fast[which(Icand2fast==0)] <-  max(Icand2fast)}
  
  which.min(Icand1fast)
  which.min(Icand2fast)
  
  ### Fast update; Equation 6.6. in Surrogates ###
  ### ALC; How much can be improved. Equation 6.6. in Surrogates ###
  alcfast <- c(Icand1fast[which.min(Icand1fast)], Icand2fast[which.min(Icand2fast)])
  alcfast
  
  ### cost; 1, 2, 3 ###
  which.min(alcfast*c(0.62, (0.62+1.03)))
  alcfast*c(0.62, (0.62+1.03))
  
  
  chosen <- rbind(chosen, c(which.min(alcfast*c(0.62, (0.62+1.03))), which.min(cbind(Icand1fast, Icand2fast)[,which.min(alcfast*c(0.62, (0.62+1.03)))])))
  Iselect <- IMSPEselect2(X.test, X.test[chosen[nrow(chosen),2],], Iselect$fit, level=chosen[nrow(chosen),1])
  
  nonlinear.cost
  nonlinear.error
  
  if(nonlinear.cost[length(nonlinear.cost)] >= 30){break}
  
}


### Save results ###
nonlinear.cost
nonlinear.error


> nonlinear.cost
[1]  0.00  1.65  3.30  4.95  6.60  8.25  9.90
[8] 11.55 13.20 14.85 16.50 18.15 19.80 21.45
[15] 23.10 24.75 25.37 
> nonlinear.error
[1] 0.42017537 0.09080948 0.08437654
[4] 0.09479839 0.09301535 0.09336047
[7] 0.10123101 0.10607513 0.10600256
[10] 0.09898270 0.09822202 0.09559912
[13] 0.09737082 0.10078014 0.10286895
[16] 0.09151262 0.08884097 

