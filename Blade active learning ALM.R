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


### RMSE ###
sqrt(sum((predy-y.test)^2))/(sqrt(sum((y.test)^2))) # closed form


### IMSPE ###
Icurrent <- mean(predsig2) # current IMSPE
Icurrent

### Add 1 points and calculate IMSPE ###
Icand1fast <- IMSPEfast2(X.test, X.test[which.max(predsig2),], fit.closed, level=1)$IMSPE
Icand2fast <- IMSPEfast2(X.test, X.test[which.max(predsig2),], fit.closed, level=2)$IMSPE

alcfast <- c(Icurrent - Icand1fast, Icurrent - Icand2fast)
alcfast

### cost; 1, 2, 3 ###
which.max(alcfast/c(0.62, (0.62+1.03)))
alcfast/c(0.62, (0.62+1.03))


chosen <- matrix(0, ncol=2)
chosen[1,1] <- which.max(alcfast/c(0.62, (0.62+1.03)))
chosen[1,2] <- which.max(predsig2)


nonlinear.cost <- 0
nonlinear.error <- sqrt(sum((predy-y.test)^2))/(sqrt(sum((y.test)^2)))

Iselect <- IMSPEselect2(X.test, X.test[chosen[nrow(chosen),2],], fit.closed, level=chosen[nrow(chosen),1]) 



#################
### Add point ###
#################
while(nonlinear.cost[length(nonlinear.cost)] < 25){ # if total cost is less than the budget
  
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
  Icand1fast <- IMSPEfast2(X.test, X.test[which.max(predsig2),], fit.closed, level=1)$IMSPE
  Icand2fast <- IMSPEfast2(X.test, X.test[which.max(predsig2),], fit.closed, level=2)$IMSPE

  alcfast <- c(Icurrent - Icand1fast, Icurrent - Icand2fast)
  alcfast
  
  which.max(alcfast/c(0.62, (0.62+1.03)))
  alcfast/c(0.62, (0.62+1.03))
  
  
  chosen <- rbind(chosen, c(which.max(alcfast/c(0.62, (0.62+1.03))), which.max(predsig2)))
  Iselect <- IMSPEselect2(X.test, X.test[chosen[nrow(chosen),2],], Iselect$fit, level=chosen[nrow(chosen),1])
  
  nonlinear.cost
  nonlinear.error
  
  if(nonlinear.cost[length(nonlinear.cost)] >= 25){break}
  
}


### Save results ###
nonlinear.cost
nonlinear.error


> nonlinear.cost
[1]  0.00  1.65  3.30  4.95  6.60  8.25  9.90
[8] 11.55 13.20 14.85 16.50 18.15 19.80 21.45
[15] 23.10 24.75 26.40
> nonlinear.error
[1] 0.42017537 0.09141011 0.07731295
[4] 0.08687695 0.07247368 0.05709499
[7] 0.05530500 0.05932019 0.05663349
[10] 0.04238410 0.03848001 0.04574467
[13] 0.04471993 0.04522004 0.03197458
[16] 0.02974547 0.03033203

