source('/Users/junoh/Downloads/code_technometrics/Internal_functions.R')
source('/Users/junoh/Downloads/code_technometrics/Sequential_algorithms.R')
library('DiceKriging')
library('rgenoud')
library('mcmc')
install.packages("devtools")
library(devtools)
install_github("cran/MuFiCokriging")
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
eps <- sqrt(.Machine$double.eps)

#-----------------------------------------------------------------------------#
#----------------#             sequential co-kriging         #----------------#
#-----------------------------------------------------------------------------#

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


#--- Initial experimental design sets
n1 <- 10; n2 <- 5

set.seed(2) 

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



#--- Multi-fidelity co-kriging model building
mymodel <- MuFicokm(
  formula = list(~1,~1), 
  MuFidesign = NestDesign, 
  response = list(y1,y2), 
  lower=eps, upper=0.1,
  nlevel = 2,
  covtype = "gauss"
)
res <- summary(mymodel)
#--- Sequential design
predictions <- predict(
  object = mymodel, 
  newdata = X.test,
  type="SK")
###
rmseco <- sqrt(sum((predictions$mean-y.test)^2))/(sqrt(sum((y.test)^2))) # Cokm
costco <- c(0)

## One point at-a-time sequential cokriging (see Section 3.1). ##
B <- 1.03/0.62		#-B: ratio of computational costs between level 1 and 2

# for(i in 1:25){
  ###
  # set.seed(9) ### 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
  # 
  # n1 <- 12; n2 <- 9
  # 
  # X1 <- maximinLHS(n1, 1)
  # X2 <- maximinLHS(n2, 1)
  # 
  # NestDesign <- NestedDesignBuild(design = list(X1,X2))
  # 
  # X1 <- NestDesign$PX
  # X2 <- ExtractNestDesign(NestDesign,2)
  # 
  # y1 <- f1(X1)
  # y2 <- f2(X2)
  # 
  # #--- Multi-fidelity co-kriging model building
  # mymodel <- MuFicokm(
  #   formula = list(~1,~1), 
  #   MuFidesign = NestDesign, 
  #   response = list(y1,y2), 
  #   # lower=eps, upper=0.1,
  #   nlevel = 2,
  #   covtype = "gauss"
  # )
  # res <- summary(mymodel)
  # #--- Sequential design
  # predictions <- predict(object = mymodel, newdata = X,test, type="SK")
  # 
  ###

for(i in 1:41){
  
  
  set.seed(1)
  niter <- 14
  
  cokm_varmax <- one_step_cokm_varmax(
    model = mymodel,
    B = B,
    xpred = X.test,
    yreal = y.test,
    # myfunctions = list(f1,f2),
    niter = niter,
    param.estim = TRUE,
    error.compute = FALSE,
    error.LOO = FALSE,
    ponderation = FALSE)
  
  ###
  rmseco <- c(rmseco, sqrt(sum((cokm_varmax$ypredseq$mean-y.test)^2))/(sqrt(sum((y.test)^2)))) # Cokm
  costco <- c(0, 0.62*cokm_varmax$CoutSave)
  rmseco
  costco
  
  
  }
  # }
  
cokm_varmax$CoutSave <- 0.62*cokm_varmax$CoutSave
###
cokm_varmax$CoutSave
rmseco
costco

costco <- c(0.00, 0.62, 1.24, 1.86, 2.48, 3.10, 3.72,
            4.34, 4.96, 5.58, 6.20, 6.82, 7.44, 8.06,
            8.68, 9.30, 9.92, 10.54, 11.16, 11.78, 12.40,
            13.02, 13.64, 14.26, 14.88, 15.50, 16.12, 16.74,
            17.36, 17.98, 18.60, 19.22, 19.84, 20.46, 21.08,
            21.70, 22.32, 22.94, 23.56, 24.18, 24.80, 25.42)

rmseco <- c(0.32067186, 0.08101806, 0.05870327, 0.06147259, 0.06420197, 0.06498355, 0.06525282, 0.06655419, 0.06557848, 0.06837311,
            0.37284240, 0.06256586, 0.06107903, 0.06632035, 0.06103831, 0.04946586, 0.05572502, 0.04933187, 0.36791050, 0.34565718,
            0.34382267, 0.44097137, 0.36732516, 0.15466829, 0.41567397, 0.33651872, 0.17770533, 0.02523661, 0.33893688, 0.12748455,
            0.34887030, 0.31978720, 0.31891140, 0.31792560, 0.09888460, 0.02193176, 0.31752500, 0.30307003, 0.32668211, 0.14002788,
            0.29572191, 0.05536005)



