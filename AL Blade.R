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
fit.closed <- closed(X1, y1, X2, y2, kernel="sqex", constant=TRUE)
predy <- predclosed(fit.closed, X.test)$mu
predsig2 <- predclosed(fit.closed, X.test)$sig2


### RMSE ###
sqrt(mean((predy-y.test)^2)) # closed form


### IMSPE ###
Icurrent <- mean(predsig2) # current IMSPE
Icurrent

### Add 1 points and calculate IMSPE ###
intgvr <- integvar(X.test, fit.closed, mc.sample=100)

Icand1fast <- intgvr$intvar1
Icand2fast <- intgvr$intvar2

which.min(Icand1fast)
which.min(Icand2fast)

### Fast update; Equation 6.6. in Surrogates ###
### ALC; How much can be improved. Equation 6.6. in Surrogates ###
alcfast <- c(Icurrent - Icand1fast[which.min(Icand1fast)], Icurrent - Icand2fast[which.min(Icand2fast)])
alcfast

### cost; 1, 2, 3 ###
which.max(alcfast/c(0.62, (0.62+1.03)))
alcfast/c(0.62, (0.62+1.03))


chosen <- matrix(0, ncol=2)
chosen[1,1] <- which.max(alcfast/c(0.62, (0.62+1.03)))
chosen[1,2] <- which.min(cbind(Icand1fast, Icand2fast)[,chosen[1,1]])


nonlinear.cost <- 0
nonlinear.error <- sqrt(mean((predy-y.test)^2))

Iselect <- IMSPEselect2(X.test[chosen[nrow(chosen),2],], fit.closed, level=chosen[nrow(chosen),1])



#################
### Add point ###
#################
while(nonlinear.cost[length(nonlinear.cost)] < 25){ # if total cost is less than the budget
  
  ### closed ###
  predy <- predclosed(Iselect$fit, X.test)$mu
  predsig2 <- predclosed(Iselect$fit, X.test)$sig2
  
  ### RMSE ###  
  nonlinear.error <- c(nonlinear.error, sqrt(mean((predy-y.test)^2))) # closed form
  if(chosen[nrow(chosen),1] == 1){
    nonlinear.cost[length(nonlinear.cost)+1] <- nonlinear.cost[length(nonlinear.cost)]+0.62
  }else if(chosen[nrow(chosen),1] == 2){
    nonlinear.cost[length(nonlinear.cost)+1] <- nonlinear.cost[length(nonlinear.cost)]+(0.62+1.03)
  }
  
  #############
  ### IMSPE ###
  #############
  Icurrent <- mean(predsig2)
  
  ### Add 1 points to the low-fidelity data ###
  intgvr <- integvar(X.test, Iselect$fit, mc.sample=100)
  
  Icand1fast <- intgvr$intvar1
  Icand2fast <- intgvr$intvar2
  
  which.min(Icand1fast)
  which.min(Icand2fast)
  
  ### Fast update; Equation 6.6. in Surrogates ###
  ### ALC; How much can be improved. Equation 6.6. in Surrogates ###
  alcfast <- c(Icurrent - Icand1fast[which.min(Icand1fast)], Icurrent - Icand2fast[which.min(Icand2fast)])
  alcfast
  
  ### cost; 1, 2, 3 ###
  which.max(alcfast/c(0.62, (0.62+1.03)))
  alcfast/c(0.62, (0.62+1.03))
  
  
  chosen <- rbind(chosen, c(which.max(alcfast/c(0.62, (0.62+1.03))), which.min(cbind(Icand1fast, Icand2fast)[,which.max(alcfast/c(0.62, (0.62+1.03)))])))
  Iselect <- IMSPEselect2(X.test[chosen[nrow(chosen),2],], Iselect$fit, level=chosen[nrow(chosen),1])
  
  nonlinear.cost
  nonlinear.error
  
  if(nonlinear.cost[length(nonlinear.cost)] >= 25){break}
  
}


### Save results ###
nonlinear.cost
nonlinear.error


> nonlinear.cost
[1]  0.00  0.62  1.24  2.89  3.51  4.13  4.75
[8]  5.37  5.99  6.61  8.26  8.88 10.53 12.18
[15] 13.83 15.48 17.13 18.78 20.43 22.08 23.73
[22] 25.38
> nonlinear.error
[1] 3.2611319 1.9233651 1.8363613 1.5997293
[5] 1.3733761 1.3515449 1.3865496 1.3740616
[9] 1.3680320 1.3767876 1.2847067 1.3014988
[13] 1.1535156 0.8852040 0.8428560 0.8294780
[17] 0.7333394 0.7886327 0.7166451 0.7190991
[21] 0.6916721 0.7031986

