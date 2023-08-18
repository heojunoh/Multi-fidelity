# install.packages("devtools")
# library(devtools)
# install_github("cran/MuFiCokriging")
library(MuFiCokriging)
library(lhs)
library(plgp)

### reproducing Section 5.3: Thermal Stress Analysis of Jet Engine Turbine Blade ###
library(matlabr)
library(randtoolbox)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)
source("/Users/junoh/Downloads/StackingDesign-Reproducibility/GP.R")
source("/Users/junoh/Downloads/StackingDesign-Reproducibility/matern.R")
source("/Users/junoh/Downloads/StackingDesign-Reproducibility/stacking_design.R")

library(laGP)
source("GP.R")
source("closed.R")
source("score.R")

cost <- NULL
d <- 2              # d: dimension of X (scalar)
n.init <- 5*d
alpha <- NULL
tt <- 2             # mesh = 0.4/(tt^l), lowest level = 0.2
Lmax <- 5
k <- NULL           # k: the parameter for the Matern kernel function (NULL means it'll be estimated by LOOCV) 
n.max <- 300        # the maximum number of sample size
log.fg <- TRUE
l <- 5

eps <- sqrt(.Machine$double.eps)

rep <- 50
result.blade.rmse <- matrix(NA, rep, 2)
colnames(result.blade.rmse) <- c("closed", "Cokriging")
result.blade.meanscore <- matrix(NA, rep, 2)
colnames(result.blade.meanscore) <- c("closed", "Cokriging")
result.blade.meancrps <- matrix(NA, rep, 2)
colnames(result.blade.meancrps) <- c("closed", "Cokriging")
result.blade.comptime <- matrix(NA, rep, 2)
colnames(result.blade.comptime) <- c("closed", "Cokriging")

n1 <- 16; n2 <- 12; n3 <- 8

### Test data ###
# n <- 100
# set.seed(1)
# X.test <- maximinLHS(n, d)
# 
# # MM <- 0.3
# d1 <- data.frame(X.test*0.5+0.25, rep(0.0125, nrow(X.test))) # scale X to [-1,1]
# write.csv(d1, "/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/temp_to_matlab.txt", row.names=F)
# write.csv(X.test, "/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/Xtest.txt", row.names=F)
# write.csv(y.test, "/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/ytest.txt", row.names=F)
# run_matlab_script("/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/SolveJetBlade.m", verbose = FALSE, desktop = FALSE, 
#                   splash = FALSE, display = FALSE, wait = TRUE, single_thread = FALSE,
#                   intern = TRUE)
# d2 <- read.table("/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/temp_to_r.txt", sep = ",")
# y.test <- d2$V4 # MM=0.3 takes 10 minutes #

for(i in 1:rep) {
  
  # i <- 1
  set.seed(i)
  
  ### Generate Input ###
  X1 <- maximinLHS(n1, d)
  X2 <- maximinLHS(n2, d)
  X3 <- maximinLHS(n3, d)
  
  NestDesign <- NestedDesignBuild(design = list(X1,X2,X3))
  
  X1 <- NestDesign$PX
  X2 <- ExtractNestDesign(NestDesign,2)
  X3 <- ExtractNestDesign(NestDesign,3)
  
  ### Y1 ###
  d1 <- data.frame(X1*0.5+0.25, rep(0.05, nrow(X1))) # scale X to [-1,1]
  write.csv(X1, "/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/temp_to_X.txt", row.names=F)
  write.csv(d1, "/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/temp_to_matlab.txt", row.names=F)
  run_matlab_script("/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/SolveJetBlade.m", verbose = FALSE, desktop = FALSE, 
                    splash = FALSE, display = FALSE, wait = TRUE, single_thread = FALSE,
                    intern = TRUE)
  d2 <- read.table("/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/temp_to_r.txt", sep = ",")
  y1 <- d2$V4
  y1
  
  ### Y2 ###
  d1 <- data.frame(X2*0.5+0.25, rep(0.025, nrow(X2))) # scale X to [-1,1]
  write.csv(X2, "/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/temp_to_X.txt", row.names=F)
  write.csv(d1, "/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/temp_to_matlab.txt", row.names=F)
  run_matlab_script("/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/SolveJetBlade.m", verbose = FALSE, desktop = FALSE, 
                    splash = FALSE, display = FALSE, wait = TRUE, single_thread = FALSE,
                    intern = TRUE)
  d2 <- read.table("/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/temp_to_r.txt", sep = ",")
  y2 <- d2$V4
  y2
  
  ### Y3 ###
  d1 <- data.frame(X3*0.5+0.25, rep(0.0125, nrow(X3))) # scale X to [-1,1]
  write.csv(X3, "/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/temp_to_X.txt", row.names=F)
  write.csv(d1, "/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/temp_to_matlab.txt", row.names=F)
  run_matlab_script("/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/SolveJetBlade.m", verbose = FALSE, desktop = FALSE, 
                    splash = FALSE, display = FALSE, wait = TRUE, single_thread = FALSE,
                    intern = TRUE)
  d2 <- read.table("/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/temp_to_r.txt", sep = ",")
  y3 <- d2$V4
  y3
  
  
  ### closed ###
  tic.closed <- proc.time()[3]
  fit.closed <- closed2(X1, y1, X2, y2, X3, y3, constant=TRUE)
  predy <- predclosed2(fit.closed, X.test)$mu
  predsig2 <- predclosed2(fit.closed, X.test)$sig2
  toc.closed <- proc.time()[3]
  
  
  ### prediction of original GP with single fidelity ###
  fit.GP3 <- GP(X3, y3, constant=TRUE)
  pred3 <- pred.GP(fit.GP3, X.test)
  
  ### Cokriging ###
  tic.cokm <- proc.time()[3]
  fit.muficokm <- MuFicokm(formula = list(~1,~1,~1), MuFidesign = NestDesign, covtype="gauss",
                           lower=eps, #upper=0.1,
                           coef.trend = list(0,c(0,0),c(0,0)), response = list(y1,y2,y3), nlevel = 3)
  pred.muficokm <- predict(fit.muficokm, X.test, "SK")
  toc.cokm <- proc.time()[3]
  
  ### RMSE ###
  result.blade.rmse[i,1] <- sqrt(mean((predy-y.test)^2)) # closed form
  result.blade.rmse[i,2] <- sqrt(mean((pred.muficokm$mean-y.test)^2)) # Cokm
  
  result.blade.meanscore[i,1] <- mean(score(y.test, predy, predsig2)) # closed form
  result.blade.meanscore[i,2] <- mean(score(y.test, pred.muficokm$mean, pred.muficokm$sig2)) # Cokriging
  
  result.blade.meancrps[i,1] <- mean(crps(y.test, predy, predsig2)) # closed form
  result.blade.meancrps[i,2] <- mean(crps(y.test, pred.muficokm$mean, pred.muficokm$sig2)) # Cokriging
  
  result.blade.comptime[i,1] <- toc.closed - tic.closed
  result.blade.comptime[i,2] <- toc.cokm - tic.cokm
  
}

par(mfrow=c(1,1))
#RMSE comparison#
apply(result.blade.rmse, 2, mean)
table(apply(result.blade.rmse, 1, which.min))
boxplot(result.blade.rmse)


