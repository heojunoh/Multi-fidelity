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

crps <- function(x, mu, sig2){ # The smaller, the better (0 to infinity)
  if(any(sig2==0)) sig2[sig2==0] <- eps
  -sqrt(sig2)*(1/sqrt(pi)-2*dnorm((x-mu)/sqrt(sig2))-(x-mu)/sqrt(sig2)*(2*pnorm((x-mu)/sqrt(sig2))-1))
}

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

rep <- 100
result.blade.rmse <- matrix(NA, rep, 2)
colnames(result.blade.rmse) <- c("RNAmf", "Cokriging")
result.blade.meancrps <- matrix(NA, rep, 2)
colnames(result.blade.meancrps) <- c("RNAmf", "Cokriging")
result.blade.comptime <- matrix(NA, rep, 2)
colnames(result.blade.comptime) <- c("RNAmf", "Cokriging")

n1 <- 30; n2 <- 20; n3 <- 10

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

# for(i in 1:rep) {
  
  i <- 100
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
  fit.closed <- RNAmf2(X1, y1, X2, y2, X3, y3, kernel="sqex", constant=TRUE)
  predy <- predRNAmf2(fit.closed, X.test)$mu
  predsig2 <- predRNAmf2(fit.closed, X.test)$sig2
  toc.closed <- proc.time()[3]
  
  
  ### Cokriging ###
  tic.cokm <- proc.time()[3]
  fit.muficokm <- MuFicokm(formula = list(~1,~1,~1), MuFidesign = NestDesign, covtype="gauss",
                           lower=eps, upper=1,
                           # coef.trend = list(0,c(0,0),c(0,0)), 
                           response = list(y1,y2,y3), nlevel = 3)
  pred.muficokm <- predict(fit.muficokm, X.test, "SK")
  toc.cokm <- proc.time()[3]
  
  ### Cokriging ###
  tic.cokmc <- proc.time()[3]
  fit.muficokmc <- MuFicokm(formula = list(~1,~1,~1), MuFidesign = NestDesign, covtype="gauss",
                           lower=eps, upper=1,
                           coef.trend = list(0,c(0,0),c(0,0)), 
                           response = list(y1,y2,y3), nlevel = 3)
  pred.muficokmc <- predict(fit.muficokmc, X.test, "SK")
  toc.cokmc <- proc.time()[3]
  
  ### RMSE ###
  result.blade.rmse[i,1] <- sqrt(mean((predy-y.test)^2)) # closed form
  result.blade.rmse[i,2] <- sqrt(mean((pred.muficokm$mean-y.test)^2)) # Cokm
  
  result.blade.meancrps[i,1] <- mean(crps(y.test, predy, predsig2)) # closed form
  result.blade.meancrps[i,2] <- mean(crps(y.test, pred.muficokm$mean, pred.muficokm$sig2)) # Cokriging
  
  result.blade.comptime[i,1] <- toc.closed - tic.closed
  result.blade.comptime[i,2] <- toc.cokm - tic.cokm
  
# }


sqrt(mean((predy-y.test)^2))
mean(crps(y.test, predy, predsig2))
toc.closed - tic.closed


sqrt(mean((pred.muficokm$mean-y.test)^2))
mean(crps(y.test, pred.muficokm$mean, pred.muficokm$sig2))
toc.cokm - tic.cokm


sqrt(mean((pred.muficokmc$mean-y.test)^2))
mean(crps(y.test, pred.muficokmc$mean, pred.muficokmc$sig2))
toc.cokmc - tic.cokmc


### RNAmf ###

c(0.8847639,0.7818165,1.059549,2.552126,1.059758,0.7242456,0.8897539,1.022244,1.410776,1.664944,
  1.27399,1.44195,1.582749,1.336593,2.473515,1.382957,1.061654,1.233598,1.7289,1.178341,
  2.250495,1.437135,1.235228,0.6474449,1.275611,1.269578,1.109555,1.213696,2.475595,1.254129,
  1.230133,1.743223,2.047466,1.020915,1.988568,1.24765,0.9160801,1.150359,0.9684565,1.464755,
  0.7108508,1.541299,1.982033,1.204411,1.704594,1.968563,1.516999,1.001181,2.170426,1.003866,
  1.355939,2.619515,1.345387,1.917427,1.331748,2.209935,2.035158,1.642062,1.966994,1.690548,
  1.630144,1.060738,0.8853699,1.251579,1.721028,2.088621,1.004609,1.593533,1.371455,1.421321,
  1.631248,1.076297,1.341311,1.534833,1.125198,0.8054764,0.9808146,0.9284471,1.239097,1.500164,
  1.657319,1.634706,1.214069,0.9585306,1.167976,1.812102,1.353286,1.117416,1.237789,0.7496117,
  0.8603298,0.9914729,0.864677,1.171613,1.859671,2.221836,0.8778934,1.428683,1.457011,1.87554
)


c(0.4727261,0.3311073,0.7272317,1.113689,0.5366611,0.7250871,0.4100703,0.4773007,0.7260806,0.6332435,
  0.7447595,0.7670454,0.8184033,0.7727856,0.9673654,0.7355679,0.4890373,0.5592317,0.7459395,0.6142787,
  1.021708,0.6449517,0.6458511,0.346264,0.5599838,0.644014,0.4408466,0.5242648,0.9062168,0.7470027,
  0.5261465,0.6672963,0.7802475,0.4802666,0.7880035,0.5639716,0.4177389,0.6001778,0.4605291,0.7953163,
  0.3293173,0.6321073,0.9546654,0.5966756,0.8115007,0.7424217,0.7590241,0.6433804,1.340091,0.5749644,
  0.5453783,0.9358727,0.5946177,0.9622639,0.6712402,0.880429,0.6373563,0.5950054,0.7203377,0.6043942,
  0.7923953,0.5115618,0.4285083,0.7308884,0.9616046,1.055698,0.5311776,0.621903,0.6813285,0.6791349,
  0.8062874,0.5541856,0.6462557,0.7059647,0.8077265,0.5126866,0.5936308,0.4885819,0.6175156,0.7505712,
  0.7805226,0.6597296,0.6107067,0.4041231,0.4775676,0.7294296,0.7370933,0.4584988,0.6029775,0.3767766,
  0.5850501,0.6560437,0.4397538,0.5944593,0.9800515,0.8848724,0.4844794,0.6371504,0.7813648,0.8035827
)


c(0.924,0.616,0.949,0.776,1.074,0.955,1.003,1.015,1.064,1.26,
  1.07,1.091,1.022,1.144,0.68,0.581,0.971,0.964,0.995,1.185,
  0.679,1.275,0.976,1.055,1.105,0.714,1.419,1.202,0.697,1.041,
  1.113,1.033,1.139,1.07,1.023,1.205,1.037,1.057,1.052,1.039,
  1.001,0.706,1.112,1.03,1.014,0.957,1.062,0.994,1.136,0.945,
  1.098,0.753,1.081,1.029,1.073,1.099,1.007,1.016,1.167,0.988,
  1.025,0.972,1.064,0.263,0.939,0.719,0.664,1.236,1.073,1.047,
  0.549,1.038,1.011,1.036,1.1,1.052,1.013,1.008,0.711,1.214,
  1.029,1.045,0.984,0.982,1.059,1.052,1.132,0.951,0.964,0.976,
  1.096,1.093,1.004,0.968,0.969,0.744,1.134,1.061,1.034,1.136
)


### CoKriging functional mean ###

c(0.8633038,0.8261792,0.650414,1.206401,0.9822952,0.798501,0.6192951,0.9076432,0.7678808,0.7815889,
  0.7812057,0.9171927,0.9180921,0.7833132,1.289448,0.9707832,0.9271943,1.188823,0.5718617,0.7578455,
  1.150132,0.813186,1.073261,1.097447,0.9616484,0.9113635,0.7822783,0.8763907,0.7488283,0.8741257,
  0.5976706,0.6425069,0.9507278,1.08207,0.6347854,0.6492633,0.863336,0.912453,0.6164437,0.7731585,
  0.6487167,0.4761606,0.8702708,0.7367809,1.325684,0.8513253,1.653201,0.832677,1.468859,0.8481921,
  0.7067045,1.893825,0.6431889,1.054912,0.9031101,0.9917614,0.7768668,0.7790027,1.097563,0.8807483,
  0.8492026,1.265218,0.8423742,1.265212,1.627464,1.920611,0.7002125,1.668208,1.498427,1.155116,
  0.9053661,0.8403947,0.9096109,1.04055,1.018651,0.7238234,0.6502924,0.8625495,0.9480625,0.6364141,
  1.025592,0.6350025,1.111441,0.9017527,0.5390661,1.026303,0.8356111,0.8798574,0.8378873,0.8241218,
  0.6796647,0.8383349,0.7949806,0.6853636,0.8479397,0.9113289,0.8009649,0.9629174,0.9037593,1.561739
)


c(0.4339874,0.3651481,0.2928728,0.511577,0.4576677,0.3997547,0.3146971,0.5029579,0.3348897,0.3701597,
  0.4056562,0.4094138,0.511578,0.3906344,0.4625985,0.447145,0.3926555,0.5582494,0.2992601,0.3418951,
  0.5358509,0.4085507,0.4378003,0.3786567,0.4762616,0.4554289,0.4026424,0.3717816,0.3328952,0.3859185,
  0.3014939,0.3534377,0.3555174,0.4918688,0.2825306,0.3460937,0.4144759,0.4403043,0.3130134,0.3723139,
  0.3354318,0.2498638,0.4249529,0.361634,0.5367981,0.374155,0.7627041,0.3764195,0.7511415,0.3768397,
  0.3813373,1.050451,0.2700536,0.431369,0.3934708,0.4346927,0.3927212,0.4122908,0.5073769,0.3821448,
  0.3647113,0.5738822,0.4133764,0.7010996,0.9009932,0.9337278,0.3276596,0.5682021,0.5241694,0.5100793,
  0.4223921,0.4177108,0.4087069,0.4283387,0.4330835,0.3759831,0.333504,0.4232065,0.4782989,0.327058,
  0.4698988,0.2740712,0.4581808,0.4180861,0.2339979,0.5020739,0.4177907,0.3643035,0.4665285,0.4140048,
  0.3361008,0.4053165,0.43672,0.3211744,0.4028885,0.4502088,0.4494417,0.4997963,0.4122696,0.6685496
)


c(0.109,0.077,0.083,0.093,0.095,0.073,0.08,0.083,0.073,0.086,
  0.1,0.08,0.079,0.155,0.084,0.083,0.072,0.079,0.087,0.077,
  0.084,0.08,0.078,0.08,0.077,0.098,0.089,0.078,0.087,0.081,
  0.094,0.075,0.082,0.094,0.098,0.105,0.103,0.123,0.074,0.091,
  0.087,0.077,0.075,0.075,0.082,0.083,0.092,0.076,0.088,0.073,
  0.084,0.099,0.078,0.087,0.075,0.078,0.099,0.078,0.088,0.088,
  0.101,0.09,0.095,0.069,0.082,0.082,0.078,0.08,0.08,0.093,
  0.076,0.164,0.075,0.078,0.098,0.072,0.072,0.089,0.079,0.076,
  0.08,0.073,0.073,0.095,0.084,0.078,0.096,0.103,0.073,0.076,
  0.1,0.078,0.076,0.106,0.083,0.075,0.072,0.091,0.091,0.096
)


### CoKriging constant mean ###

c(1.741914,1.848022,2.325056,3.449001,2.259547,2.540847,3.583533,2.325992,2.782767,2.910051,
  3.626293,3.264911,2.153301,2.348823,2.662095,2.463939,1.938478,2.081244,1.820247,3.941962,
  3.753655,2.16167,1.7454,2.232516,3.395771,1.593151,1.863782,5.369144,2.216523,2.917966,
  1.726851,3.458678,2.47178,1.887852,2.477598,2.560552,2.509144,2.124056,1.885385,2.122871,
  2.081,2.678941,3.335142,1.849865,3.950518,3.423523,2.281278,3.79985,2.961052,2.192837,
  2.062142,2.331321,3.061982,2.601393,2.554739,2.508237,2.582427,2.456111,2.808877,3.26927,
  1.951511,2.436778,4.481934,2.261543,1.463976,1.943478,1.995274,2.208513,2.068537,1.62307,
  1.613802,1.863823,5.079342,1.835687,5.109532,4.444028,2.111074,2.048751,1.687471,2.532323,
  2.884717,2.318962,2.018802,1.728134,2.007855,3.046823,2.156287,2.351164,1.563071,1.721056,
  1.463296,2.840496,2.012369,2.123552,1.766857,3.40627,1.870899,3.190091,1.816971,2.434534
)


c(0.9450146,1.113146,1.17276,1.747752,1.119361,1.437096,1.923035,1.385822,1.223689,1.678157,
  1.879902,1.385041,1.216437,1.464511,1.369393,1.390308,1.040633,1.137484,0.9285124,1.867006,
  1.953653,1.163544,0.904444,1.402958,2.014963,0.8545805,0.9679488,3.039253,1.150668,1.411885,
  0.9055712,1.875515,1.322309,1.019236,1.274587,1.244258,1.377255,1.222493,0.9995321,1.044145,
  1.14133,1.45266,2.012585,1.20251,2.142988,1.67962,1.240877,1.837226,1.496853,1.169286,
  1.030658,1.303268,1.486274,1.385251,1.641541,1.401308,1.39619,1.249895,1.419121,1.395592,
  1.113226,1.381917,2.162022,1.488748,0.8372623,1.085758,1.130079,1.029485,1.120628,0.9194718,
  0.903592,1.100474,2.958282,1.117789,2.999873,2.491565,1.188309,1.154667,0.8557967,1.592126,
  1.586307,1.080808,1.14575,0.9200112,1.167994,1.412776,1.256306,1.33478,0.8014271,1.009874,
  0.8465369,1.759105,1.05228,1.179777,0.9570332,1.543012,0.9759724,1.422411,1.089633,1.167494
)


c(0.079,0.067,0.062,0.071,0.092,0.061,0.061,0.063,0.107,0.068,
  0.083,0.064,0.099,0.073,0.07,0.06,0.061,0.065,0.065,0.073,
  0.068,0.072,0.065,0.076,0.071,0.077,0.087,0.096,0.076,0.072,
  0.068,0.067,0.071,0.066,0.088,0.07,0.097,0.069,0.069,0.066,
  0.071,0.076,0.071,0.066,0.066,0.065,0.063,0.068,0.072,0.077,
  0.066,0.063,0.067,0.077,0.088,0.075,0.068,0.064,0.073,0.065,
  0.083,0.065,0.064,0.06,0.056,0.058,0.068,0.062,0.072,0.065,
  0.079,0.062,0.065,0.068,0.075,0.073,0.075,0.062,0.067,0.072,
  0.081,0.062,0.065,0.066,0.065,0.081,0.072,0.076,0.075,0.061,
  0.076,0.062,0.058,0.066,0.064,0.065,0.07,0.068,0.076,0.067
)




# par(mfrow=c(1,1))
# #RMSE comparison#
# apply(result.blade.rmse, 2, mean)
# table(apply(result.blade.rmse, 1, which.min))
# boxplot(result.blade.rmse)

