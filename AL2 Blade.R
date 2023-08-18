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
intgvr <- integvar(matrix(X.test[which.max(predsig2),], nrow=1), fit.closed, mc.sample=100)

Icand1fast <- intgvr$intvar1
Icand2fast <- intgvr$intvar2

alcfast <- c(Icurrent - Icand1fast, Icurrent - Icand2fast)
alcfast

### cost; 1, 2, 3 ###
which.max(alcfast/c(0.62, (0.62+1.03)))
alcfast/c(0.62, (0.62+1.03))


chosen <- matrix(0, ncol=2)
chosen[1,1] <- which.max(alcfast/c(0.62, (0.62+1.03)))
chosen[1,2] <- which.max(predsig2)


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
  intgvr <- integvar(matrix(X.test[which.max(predsig2),], nrow=1), Iselect$fit, mc.sample=100)
  
  Icand1fast <- intgvr$intvar1
  Icand2fast <- intgvr$intvar2
  
  
  alcfast <- c(Icurrent - Icand1fast, Icurrent - Icand2fast)
  alcfast
  
  which.max(alcfast/c(0.62, (0.62+1.03)))
  alcfast/c(0.62, (0.62+1.03))
  
  
  chosen <- rbind(chosen, c(which.max(alcfast/c(0.62, (0.62+1.03))), which.max(predsig2)))
  Iselect <- IMSPEselect2(X.test[chosen[nrow(chosen),2],], Iselect$fit, level=chosen[nrow(chosen),1])
  
  nonlinear.cost
  nonlinear.error
  
  if(nonlinear.cost[length(nonlinear.cost)] >= 25){break}
  
}


### Save results ###
nonlinear.cost
nonlinear.error


> nonlinear.cost
[1]  0.00  0.62  1.24  2.89  3.51  4.13  5.78
[8]  6.40  7.02  7.64  8.26  8.88  9.50 11.15
[15] 11.77 12.39 13.01 13.63 14.25 14.87 15.49
[22] 16.11 16.73 17.35 17.97 18.59 19.21 19.83
[29] 20.45 21.07 21.69 22.31 22.93 23.55 24.17
[36] 24.79 25.41
> nonlinear.error
[1] 3.2611319 1.9315585 1.6128726 1.4028079
[5] 1.3854075 1.4223574 1.3481041 1.3464402
[9] 1.3205793 1.1906808 1.0899642 0.9940933
[13] 1.0188526 1.1885785 0.9822382 1.0201724
[17] 1.0693308 1.0616464 0.8565424 0.6788938
[21] 0.5354022 0.4947483 0.5140573 0.4683984
[25] 0.4608091 0.4805101 0.4781361 0.5405246
[29] 0.5066524 0.3588916 0.3694267 0.3407210
[33] 0.3538892 0.3643914 0.3408677 0.3353673
[37] 0.3198383

